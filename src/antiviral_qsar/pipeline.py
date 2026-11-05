from __future__ import annotations

import csv
import json
import re
import time
from collections import Counter, defaultdict

from rich.console import Console

from .chemistry import validate_measurement
from .config import Config, PROMPT_VERSION, RULES_VERSION, digest
from .documents import jats_chunks, zip_chunks, relevant_chunk, attachment_chunks, Chunk
from .model import LocalModel, grounded, normalized_text
from .network import Network
from .progress import RunProgress
from .schema import AssayContext, Measurement, Virus
from .sources import ChEMBL, PubChem, EuropePMC, paper_identity, pubchem_value, pubchem_unit, search_terms
from .storage import Store, now

console = Console()


def known_context(description, cell="", organism=""):
    """Only clear explicit cases bypass the model; no confidence-score cutoff for organism assays."""
    text = description.lower()
    biochemical = bool(re.search(r"surface plasmon|purified (?:enzyme|protein)|cell.free|enzyme.only|binding affinity|dissociation rate", text))
    infected = bool(re.search(r"infected|viral replication|plaque reduction|cytopathic|antiviral activity", text))
    if biochemical and not re.search(r"infected.*cells|cells.*infected", text):
        return AssayContext(kind="biochemical", evidence=description, supported=True)
    if not cell:
        names = re.findall(r"\b(?:Vero(?:\s+E6)?|Huh[ -]?7(?:\.5)?|BHK[ -]?21|A549|HEK[ -]?293(?:T)?|"
                           r"C6/36|MDCK|HepG2|MT[ -]?4|CEM|RD|HeLa)\b", description, re.I)
        names = list(dict.fromkeys(names))
        if len(names) == 1:
            cell = names[0]
    if cell and re.search(r"cytotoxic|cell viability", text) and not infected:
        return AssayContext(kind="cytotoxicity", cell_line=cell, evidence=description, supported=True)
    if cell and infected and not biochemical:
        return AssayContext(kind="cellular_antiviral", cell_line=cell, virus=organism or None,
                            evidence=description, supported=True)
    return None


def relevant_virus(text, virus):
    text = normalized_text(text)
    return any(normalized_text(term) in text for term in search_terms(virus))


class Pipeline:
    def __init__(self, config: Config, store: Store, *, use_model=True):
        self.config, self.store = config, store
        self.net = Network(config)
        self.chembl, self.pubchem, self.epmc = ChEMBL(self.net), PubChem(self.net), EuropePMC(self.net)
        self.model = LocalModel(config, store)
        self.use_model = use_model
        self.progress = RunProgress(console)

    def context(self, description, cell="", organism=""):
        result = known_context(description, cell or "", organism or "")
        if result is not None:
            return result
        if self.use_model:
            result = self.model.classify(description, cell, organism)
            result.classified_by = "llm"
            return result
        return AssayContext(cell_line=cell or None, virus=organism or None, evidence=description)

    def discover(self, rid, virus: Virus, *, assay_limit=None, paper_search_limit=None):
        inventory = self.store.run(rid)["inventory"]
        inventory.setdefault("sources", {})
        inventory["search_terms"] = search_terms(virus)
        inventory["scope"] = {"sources":["chembl", "pubchem", "europepmc"],
                              "not_connected":["AntiviralDB", "DrugVirus.info", "Stanford CoV-RDB", "BindingDB direct"],
                              "cellular_only":True, "small_molecules_only":True}
        for source in ("chembl", "pubchem", "europepmc"):
            if inventory["sources"].get(source, {}).get("status") == "available":
                continue
            console.print(f"[cyan]Searching: {source}[/cyan]")
            try:
                if source == "chembl":
                    assays, status = self.chembl.discover(virus)
                    selected = assays if assay_limit is None else assays[:assay_limit]
                    self.store.add_jobs(rid, source, "assay", ((a["assay_chembl_id"], a) for a in selected))
                    status["scheduled"] = len(selected)
                elif source == "pubchem":
                    ids, status = self.pubchem.discover(virus)
                    selected = ids if assay_limit is None else ids[:assay_limit]
                    self.store.add_jobs(rid, source, "assay", ((aid, {"aid":aid}) for aid in selected))
                    status["scheduled"] = len(selected)
                else:
                    papers, status = self.epmc.discover(virus)
                    selected = papers if paper_search_limit is None else papers[:paper_search_limit]
                    self.store.add_jobs(rid, source, "paper", ((paper_identity(p), p) for p in selected))
                    for paper in selected:
                        self.store.reuse_paper(rid, paper_identity(paper), paper, self.config.extractor_signature)
                    status["open_fulltext"] = sum(p.get("isOpenAccess") == "Y" and bool(p.get("pmcid")) for p in selected)
                    status["scheduled_candidates"] = len(selected)
                    status["previously_processed"] = sum(self.store.previously_processed(source, paper_identity(p)) for p in selected)
                inventory["sources"][source] = status
            except Exception as exc:
                inventory["sources"][source] = {"status":"unavailable", "error":str(exc)[:500]}
                self.store.event(rid, f"{source}: {exc}", "error")
                console.print(f"[yellow]{source}: unavailable ({type(exc).__name__})[/yellow]")
            self.store.update_run(rid, inventory=inventory)
        inventory["discovery_finished_at"] = now()
        inventory["sampling"] = {"assay_limit_per_source":assay_limit, "paper_search_limit":paper_search_limit}
        self.store.update_run(rid, inventory=inventory, status="awaiting_selection")
        return inventory

    def chembl_assay(self, assay, virus):
        aid = assay["assay_chembl_id"]
        rows = list(self.chembl.activities(aid))
        if not rows:
            return []
        description = assay.get("description", "")
        organism = assay.get("assay_organism") or ""
        organism_is_virus = relevant_virus(organism, virus)
        context = self.context(description, assay.get("assay_cell_type"), organism if organism_is_virus else "")
        if context.kind in ("biochemical", "other"):
            return []
        doc_id = assay.get("document_chembl_id", "")
        document = self.chembl.document(doc_id)
        origin = f"chembl:{aid}"
        if assay.get("src_id") == 7 and assay.get("src_assay_id"):
            origin = f"pubchem:{assay['src_assay_id']}"
        observations = []
        try:
            self.chembl.prefetch_molecules([a.get('molecule_chembl_id', '') for a, _ in rows])
        except Exception:
            pass  # Individual lookups below retain explicit per-assay failures.
        for activity, snapshot in rows:
            compound_id = activity.get("molecule_chembl_id", "")
            molecule = self.chembl.molecule(compound_id)
            endpoint = activity.get("standard_type", "")
            evidence = json.dumps({k:activity.get(k) for k in ("activity_id", "molecule_chembl_id", "standard_type",
                                  "standard_value", "standard_units", "standard_relation", "value", "units")}, ensure_ascii=False)
            original_value = activity.get("value")
            original_units = activity.get("units")
            original_endpoint = activity.get("type")
            # Preserve deposited values when they represent the same endpoint; standard fields remain in evidence.
            if original_value is None or not original_units or original_endpoint != endpoint:
                original_value, original_units = activity.get("standard_value"), activity.get("standard_units")
                relation = activity.get("standard_relation") or "="
            else:
                relation = activity.get("relation") or activity.get("standard_relation") or "="
            m = Measurement(observation_id=digest(f"chembl:{activity['activity_id']}"), source="chembl",
                source_record_id=str(activity["activity_id"]), source_url=f"https://www.ebi.ac.uk/chembl/explore/assay/{aid}",
                source_snapshot=snapshot, origin_assay_id=origin, document_id=doc_id, doi=document.get("doi") or "",
                compound=molecule.get("pref_name") or compound_id, compound_id=compound_id,
                smiles_original=(molecule.get("molecule_structures") or {}).get("canonical_smiles") or "",
                virus=(organism if organism_is_virus else virus.name) if endpoint != "CC50" else "",
                virus_tax_id=str(assay.get("assay_tax_id") if organism_is_virus and assay.get("assay_tax_id")
                                 else virus.tax_id) if endpoint != "CC50" else "",
                strain=assay.get("assay_strain") or context.strain or "", cell_line=context.cell_line or "",
                endpoint=endpoint, value_original=str(original_value if original_value is not None else ""),
                units_original=original_units or "", relation=relation,
                duration=context.duration or "", method=context.method or "", assay_description=description,
                evidence=evidence, evidence_locator=f"activity_id:{activity['activity_id']}", context_evidence=context.evidence,
                evidence_grounded=True, cell_context_supported=context.supported,
                context_method=context.classified_by,
                extractor_signature=self.config.extractor_signature)
            if endpoint != "CC50" and not relevant_virus(description+" "+(assay.get("assay_organism") or ""), virus):
                m.review_reasons.append("virus_scope_not_explicit")
            if endpoint != "CC50" and context.kind != "cellular_antiviral":
                m.review_reasons.append("antiviral_assay_context_unclear")
            if molecule.get("molecule_type") not in ("Small molecule", None):
                m.review_reasons.append("outside_small_molecule_scope")
            if activity.get("data_validity_comment"):
                m.warnings.append("source_validity:"+activity["data_validity_comment"])
                if activity["data_validity_comment"] != "Manually validated":
                    m.review_reasons.append("source_validity_flag")
            if activity.get("potential_duplicate"):
                m.warnings.append("source_potential_duplicate")
            observations.append(validate_measurement(m))
        return observations

    def pubchem_assay(self, aid, virus):
        description, desc_snapshot = self.pubchem.description(aid)
        definitions = {r["tid"]:r for r in description.get("results", [])}
        endpoint_defs = {tid:info for tid,info in definitions.items()
                         if re.search(r"\b(?:IC|EC|CC)\s*50\b", info.get("name", ""), re.I)}
        if not endpoint_defs and not any("standard value" in d.get("name", "").lower() for d in definitions.values()):
            return []
        text = description.get("name", "") + "\n" + "\n".join(description.get("description", []) + description.get("protocol", []))
        context = self.context(text)
        if context.kind in ("biochemical", "other"):
            return []
        source = description.get("aid_source", description.get("source", {})).get("db", {})
        origin = f"pubchem:{aid}"
        if source.get("name", "").lower() == "chembl":
            sid = source.get("source_id", {}).get("str", "")
            match = re.search(r"CHEMBL\d+", sid)
            if match:
                origin = "chembl:" + match.group()
        document_ids = []
        for entry in description.get("xref", []):
            entry = entry.get("xref", entry)
            if "pmid" in entry:
                document_ids.append(str(entry["pmid"]))
        observations = []
        rows = list(self.pubchem.results(aid))
        try:
            self.pubchem.prefetch_substances([r['sid'] for r, _ in rows])
        except Exception:
            pass  # Fallback preserves an explicit lookup error for each unresolved structure.
        for row, snapshot in rows:
            values = {field["tid"]:pubchem_value(field) for field in row.get("data", [])}
            candidates = [(tid, info) for tid,info in endpoint_defs.items() if values.get(tid)]
            if not candidates:
                standard_type = next((values.get(tid, "") for tid,d in definitions.items()
                                      if d.get("name", "").lower() == "standard type"), "")
                if standard_type in ("IC50", "EC50", "CC50"):
                    candidates = [(tid, {**info, "name":standard_type}) for tid,info in definitions.items()
                                  if "standard value" in info.get("name", "").lower() and values.get(tid)]
            if not candidates:
                continue
            sid = row["sid"]
            try:
                compound, compound_error = self.pubchem.sid_compound(sid)
            except Exception as exc:
                compound, compound_error = {}, f"compound_lookup_failed:{type(exc).__name__}"
            for tid, info in candidates:
                endpoint = re.search(r"(?:IC|EC|CC)\s*50", info["name"], re.I).group().replace(" ", "").upper()
                related = next((values.get(t, "=") for t,d in definitions.items()
                                if "relation" in d.get("name", "").lower()), "=")
                evidence = json.dumps({"aid":aid, "sid":sid, "tid":tid, "definition":info, "value":values[tid]}, ensure_ascii=False)
                m = Measurement(observation_id=digest(f"pubchem:{aid}:{sid}:{tid}"), source="pubchem",
                    source_record_id=f"{aid}:{sid}:{tid}", source_url=f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{aid}",
                    source_snapshot=snapshot, origin_assay_id=origin,
                    document_id=";".join("PMID:"+p for p in document_ids),
                    compound=f"CID:{compound.get('CID')}" if compound else f"SID:{sid}",
                    compound_id=f"CID:{compound.get('CID')}" if compound else f"SID:{sid}",
                    smiles_original=compound.get("SMILES") or compound.get("IsomericSMILES") or "",
                    virus=(context.virus or virus.name) if endpoint != "CC50" else "",
                    virus_tax_id=virus.tax_id if endpoint != "CC50" else "", cell_line=context.cell_line or "",
                    strain=context.strain or "", endpoint=endpoint, value_original=values[tid],
                    units_original=info.get("sunit") or pubchem_unit(info.get("unit")), relation=related,
                    duration=context.duration or "", method=context.method or "", assay_description=text,
                    evidence=evidence, evidence_locator=f"AID:{aid}/SID:{sid}/TID:{tid}", context_evidence=context.evidence,
                    evidence_grounded=True, cell_context_supported=context.supported,
                    context_method=context.classified_by,
                    extractor_signature=self.config.extractor_signature)
                if compound_error:
                    m.review_reasons.append(compound_error)
                if description.get("is_panel") or info.get("panel_info"):
                    m.review_reasons.append("panel_readout_requires_context_review")
                if info.get("transform", 1) != 1 or info.get("stransform"):
                    m.review_reasons.append("transformed_readout_requires_review")
                    m.value_nm = None
                if re.search(r"error|deviation|confidence|lower|upper|log|ratio", info["name"], re.I):
                    m.review_reasons.append("derived_readout_requires_review")
                if row.get("comment"):
                    m.warnings.append("source_row_comment:"+row["comment"])
                if endpoint != "CC50" and (context.kind != "cellular_antiviral" or not relevant_virus(text, virus)):
                    m.review_reasons.append("virus_or_antiviral_context_unclear")
                if context.virus and endpoint != "CC50" and not relevant_virus(context.virus, virus):
                    m.review_reasons.append("extracted_virus_scope_requires_review")
                observations.append(validate_measurement(m))
        return observations

    def extract_chunk(self, chunk, virus, depth=0):
        # Reduce a chunk instead of asking the model to silently omit measurements.
        try:
            return [(item, chunk) for item in self.model.extract(chunk.text, chunk.context, virus.name).measurements]
        except ValueError:
            lines = chunk.text.splitlines()
            if depth >= 6 or len(chunk.text) < 240:
                raise
            if len(lines) < 4:
                middle = len(chunk.text) // 2
                split = chunk.text.rfind(' ', 0, middle)
                if split < len(chunk.text) // 4:
                    split = middle
                first = Chunk(chunk.locator+':a', chunk.text[:split], chunk.context)
                second = Chunk(chunk.locator+':b', chunk.text[split:], chunk.context)
                return self.extract_chunk(first, virus, depth+1) + self.extract_chunk(second, virus, depth+1)
            middle = len(lines) // 2
            context = chunk.context[:6000]
            first = Chunk(chunk.locator+":a", "\n".join(lines[:middle]), context)
            second = Chunk(chunk.locator+":b", "\n".join(lines[middle:]), context + "\n" + "\n".join(lines[:3]))
            return self.extract_chunk(first, virus, depth+1) + self.extract_chunk(second, virus, depth+1)

    def paper(self, paper, virus, rid):
        raw, snapshot = self.epmc.fulltext(paper)
        if raw.lstrip().startswith(b'%PDF-'):
            chunks, warnings = attachment_chunks('article.pdf', raw)
            if not chunks:
                raise ValueError('main_pdf_has_no_extractable_text')
            context = (paper.get('title', '')+'\n'+'\n'.join(c.text for c in chunks[:3]))[:6000]
            for chunk in chunks:
                chunk.context = context
            for warning in warnings:
                self.store.event(rid, warning, 'warning')
            attachments = []
        else:
            chunks, attachments = jats_chunks(raw)
        combined = '\n'.join(c.text for c in chunks)
        aliases = re.findall(re.escape(virus.name) + r'\s*\(([A-Za-z][A-Za-z0-9-]{1,11})\)', combined, re.I)
        aliases += re.findall(r'\b([A-Z][A-Z0-9-]{1,11})\s*\(' + re.escape(virus.name) + r'\)', combined, re.I)
        virus = virus.model_copy(update={'synonyms':sorted(set(virus.synonyms + aliases))})
        label_names = defaultdict(set)
        for chunk in chunks:
            for name, label in re.findall(r'ROW \d+:\s*([^|()\n]{3,120})\(\s*(\d+[a-z]?)\s*\)\s*\|', chunk.text, re.I):
                label_names[label.casefold()].add(name.strip())
        errors = []
        try:
            documents = self.chembl.documents_by_reference(paper.get("doi", ""), paper.get("id", ""))
            chembl_doc = documents[0]["document_chembl_id"] if len(documents) == 1 else ""
        except Exception as exc:
            chembl_doc = ""
            self.store.event(rid, f"document_structure_lookup_unavailable:{type(exc).__name__}", "warning")
        if attachments or paper.get("hasSuppl") == "Y":
            try:
                data, supplement_snapshot = self.epmc.supplements(paper)
                if data:
                    extra, attachment_errors = zip_chunks(data, self.config.settings["collection"]["max_attachment_mb"] * 1_000_000)
                    for chunk in extra:
                        chunk.context = "Article: " + paper.get("title", "") + "\n" + (chunks[0].context if chunks else "")
                        chunk.locator = supplement_snapshot + "#" + chunk.locator
                    chunks.extend(extra)
                    for error in attachment_errors:
                        self.store.event(rid, error, "warning")
            except Exception as exc:
                error = f"supplements_unavailable:{paper.get('pmcid')}:{type(exc).__name__}"
                self.store.event(rid, error, "warning")
                errors.append(error)
        selected = [c for c in chunks if relevant_chunk(c)]
        cap = self.config.settings["collection"].get("max_article_chunks", 0)
        if cap and len(selected) > cap:
            self.store.event(rid, f"article_chunk_limit:{paper.get('pmcid')}:{len(selected)}->{cap}", "warning")
            errors.append("article_chunk_limit")
            selected = selected[:cap]
        observations = []
        for chunk in self.progress.track(selected, 'Fragments'):
            self.progress.log('fragment_started', document=paper.get('pmcid'), locator=chunk.locator)
            try:
                items = self.extract_chunk(chunk, virus)
            except Exception as exc:
                errors.append(f"{chunk.locator}:{type(exc).__name__}:{str(exc)[:180]}")
                self.progress.log('fragment_failed', locator=chunk.locator, error=str(exc))
                continue
            for item, used_chunk in items:
                evidence_ok = grounded(item.evidence, used_chunk.text)
                context_source = used_chunk.context + "\n" + used_chunk.text
                cell_ok = bool(item.cell_line and grounded(item.cell_line, context_source) and
                               grounded(item.context_evidence, context_source))
                compound, structure_error = {}, ""
                compound_id = ""
                smiles = item.smiles or ""
                if smiles and not grounded(smiles, context_source):
                    smiles, structure_error = "", "smiles_not_present_in_source"
                if not smiles:
                    try:
                        mapped = self.chembl.compound_in_document(chembl_doc, item.compound)
                        if len(mapped) == 1:
                            molecule = self.chembl.molecule(mapped[0])
                            smiles = (molecule.get("molecule_structures") or {}).get("canonical_smiles") or ""
                            compound_id = mapped[0]
                            lookup_error = "" if smiles else "unresolved_structure"
                        else:
                            name = re.sub(r'\s*\(\s*\d+[a-z]?\s*\)\s*$', '', item.compound, flags=re.I).strip()
                            explicit_names = label_names.get(name.casefold(), set())
                            if len(explicit_names) == 1:
                                name = next(iter(explicit_names))
                            compound, lookup_error = self.pubchem.resolve_name(name)
                            smiles = compound.get("SMILES") or ""
                            compound_id = f"CID:{compound['CID']}" if compound else ""
                        structure_error = structure_error or lookup_error
                    except Exception as exc:
                        structure_error = f"compound_lookup_failed:{type(exc).__name__}"
                key = [paper_identity(paper), used_chunk.locator, item.model_dump()]
                m = Measurement(observation_id=digest(key), source="europepmc", source_record_id=digest(key),
                    source_url=f"https://europepmc.org/articles/{paper.get('pmcid')}", source_snapshot=snapshot,
                    origin_assay_id=f"{paper_identity(paper)}:{item.experiment_label or used_chunk.locator}",
                    document_id=paper.get("pmcid", ""), doi=paper.get("doi", ""), compound=item.compound,
                    compound_id=compound_id, smiles_original=smiles,
                    virus=item.virus or "", virus_tax_id=virus.tax_id if item.virus and relevant_virus(item.virus, virus) else "",
                    strain=item.strain or "", cell_line=item.cell_line or "", endpoint=item.endpoint,
                    value_original=item.value, units_original=item.units, relation=item.relation,
                    duration=item.duration or "", method=item.method or "", evidence=item.evidence,
                    evidence_locator=used_chunk.locator, context_evidence=item.context_evidence,
                    assay_description=used_chunk.context, extraction_method="llm", evidence_grounded=evidence_ok,
                    context_method="llm",
                    cell_context_supported=cell_ok, extractor_signature=self.config.extractor_signature)
                if structure_error:
                    m.review_reasons.append(structure_error)
                # Token presence is necessary but not sufficient; unvalidated extraction is separately gated.
                for label, value in [("compound", item.compound), ("value", item.value.lstrip("<>= ")), ("units", item.units)]:
                    if not grounded(value, item.evidence + "\n" + used_chunk.text[:2500]):
                        m.review_reasons.append(label+"_not_grounded")
                if item.endpoint != "CC50" and not relevant_virus((item.virus or "")+" "+item.context_evidence, virus):
                    m.review_reasons.append("virus_scope_not_explicit")
                if item.endpoint == "CC50":
                    m.virus = m.virus_tax_id = ""
                observations.append(validate_measurement(m, llm_gate=self.store.quality_passed(self.config.extractor_signature)))
        return observations, errors

    def execute(self, rid, paper_budget):
        self.progress.log_path = self.config.root / 'logs' / f'{rid}.log'
        with self.progress:
            return self._execute(rid, paper_budget)

    def _execute(self, rid, paper_budget):
        run = self.store.run(rid)
        virus = Virus.model_validate(run["virus"])
        self.store.update_run(rid, status="running")
        code_snapshot = self.config.save_code_snapshot(rid)
        start = time.monotonic()
        assays = self.store.jobs(rid, "assay", unfinished=True)
        related_docs = {m.document_id for m in self.store.observations(rid)
                        if m.source == "chembl" and m.document_id}
        saved_assays = len(self.store.jobs(rid, 'assay')) - len(assays)
        for job in self.progress.track(assays, 'Assays', completed=saved_assays):
            self.progress.log('assay_started', source=job['source'], id=job['key'])
            self.store.start_job(job)
            try:
                data = (self.chembl_assay(job["payload"], virus) if job["source"] == "chembl"
                        else self.pubchem_assay(job["key"], virus))
                self.store.finish_job(job, data, "done" if data else "skipped")
                self.progress.log('assay_finished', source=job['source'], id=job['key'], records=len(data))
                if job["source"] == "chembl" and data and job["payload"].get("document_chembl_id"):
                    related_docs.add(job["payload"].get("document_chembl_id"))
            except Exception as exc:
                self.store.finish_job(job, status="failed", error=str(exc))
                self.store.event(rid, f"{job['source']}:{job['key']}:{exc}", "error")
                self.progress.log('assay_failed', source=job['source'], id=job['key'], error=str(exc))
        # Cytotoxicity often does not mention the virus; inspect companion assays from the same documents.
        if related_docs and not run["options"].get("assay_limit"):
            self.store.add_jobs(rid, "chembl", "companion_search", ((d, {"document_id":d}) for d in sorted(related_docs)))
        companion_jobs = self.store.jobs(rid, 'companion_search', unfinished=True)
        for job in self.progress.track(companion_jobs, 'Related assays'):
            self.store.start_job(job)
            try:
                for assay in self.chembl.linked_cytotoxicity([job['key']]):
                    self.store.add_job(rid, "chembl", "cytotoxicity", assay["assay_chembl_id"], assay)
                self.store.finish_job(job)
            except Exception as exc:
                self.store.finish_job(job, status='failed', error=str(exc))
                self.store.event(rid, f"companion_assays_unavailable:{job['key']}:{exc}", 'error')
        toxicity_jobs = self.store.jobs(rid, 'cytotoxicity', unfinished=True)
        for job in self.progress.track(toxicity_jobs, 'Cytotoxicity'):
            self.store.start_job(job)
            try:
                data = [m for m in self.chembl_assay(job["payload"], virus) if m.endpoint == "CC50"]
                self.store.finish_job(job, data, "done" if data else "skipped")
            except Exception as exc:
                self.store.finish_job(job, status="failed", error=str(exc))
        papers = [j for j in self.store.jobs(rid, "paper", unfinished=True)
                  if j["payload"].get("pmcid") and j["payload"].get("isOpenAccess") == "Y"]
        papers.sort(key=lambda j: (-paper_priority(j['payload'], virus), j['key']))
        if paper_budget is None and 'selected_paper_keys' in run['options']:
            keys = set(run['options']['selected_paper_keys'])
            selected = [j for j in papers if j['key'] in keys]
        else:
            budget = paper_budget if paper_budget is not None else run['options'].get('papers', '0')
            selected = papers if budget == "all" else papers[:int(budget or '0')]
            options = dict(run['options'])
            options['selected_paper_keys'] = [j['key'] for j in selected]
            self.store.update_run(rid, options=options)
        saved_papers = (len(run['options']['selected_paper_keys']) - len(selected)
                        if paper_budget is None and 'selected_paper_keys' in run['options'] else 0)
        for job in self.progress.track(selected, 'Papers', completed=saved_papers):
            self.progress.log('paper_started', id=job['key'], title=job['payload'].get('title', ''))
            self.store.start_job(job)
            try:
                data, errors = self.paper(job["payload"], virus, rid)
                self.store.finish_job(job, data, status="partial" if errors else "done", error=";".join(errors))
                self.progress.log('paper_finished', id=job['key'], records=len(data), errors=errors)
            except Exception as exc:
                self.store.finish_job(job, status="failed", error=str(exc))
                self.store.event(rid, f"paper:{job['key']}:{exc}", "error")
                self.progress.log('paper_failed', id=job['key'], error=str(exc))
        failed = [j for j in self.store.jobs(rid) if j["status"] in ("failed", "partial", "running")]
        report = self.store.run(rid)["report"]
        report["paper_selection"] = "Open texts ranked by virus in title, endpoint in abstract, primary research; ties by identifier"
        report.setdefault("sessions", []).append({"finished_at":now(), "seconds":round(time.monotonic()-start, 2),
                "code_snapshot":code_snapshot, "code_hash":self.config.code_hash,
                "llm_calls":self.model.calls, "llm_seconds":round(self.model.seconds,2),
                "prompt_tokens":self.model.tokens_in, "completion_tokens":self.model.tokens_out,
                "peak_server_rss_mb_sampled":round(self.model.peak_server_rss_mb,1),
                "network_requests":self.net.requests, "cache_hits":self.net.cache_hits})
        report["model"] = {"repo":self.config.settings["model"]["repo"], "revision":self.config.model_revision,
                           "prompt_version":PROMPT_VERSION, "rules_version":RULES_VERSION,
                           "signature":self.config.extractor_signature}
        unavailable = any(s.get("status") != "available" for s in run["inventory"].get("sources", {}).values())
        self.store.update_run(rid, status="partial" if failed or unavailable else "completed", report=report)
        return export_run(self.config, self.store, rid)


def paper_priority(paper, virus):
    terms = search_terms(virus) + [virus.name.removesuffix(" virus")]
    title = paper.get("title", "").casefold()
    abstract = paper.get("abstractText", "")
    types = " ".join(paper.get("pubTypeList", {}).get("pubType", [])).casefold()
    return (20 * any(term.casefold() in title for term in terms)
            + 10 * bool(re.search(r"\b(?:IC|EC|CC)\s*50\b", abstract, re.I))
            + 10 * bool(re.search(r"\bVero\b|Huh.?7|BHK.?21|infected cells|cell.based", abstract, re.I))
            + 5 * ("review" not in types)
            + 5 * bool(re.search(r"antiviral|inhibit|compound|molecule", title)))


def deduplicate(observations):
    groups = defaultdict(list)
    for m in observations:
        # Only merge an identical measurement tied to the same original assay, never merely same molecule.
        if m.inchikey and m.value_nm is not None and m.origin_assay_id:
            key = (m.origin_assay_id, m.inchikey, m.endpoint, m.relation, m.value_nm,
                   m.cell_line.casefold(), m.duration, m.method)
            groups[key].append(m)
    for group in groups.values():
        group.sort(key=lambda m:(len(m.review_reasons), m.source != "chembl", m.observation_id))
        # Equal repeat measurements within one source are not proven duplicates.
        if len({m.source for m in group}) != len(group):
            continue
        primary = group[0]
        for m in group:
            primary.provenance.append({"source":m.source, "record_id":m.source_record_id,
                                       "url":m.source_url, "snapshot":m.source_snapshot})
            if m is not primary:
                m.duplicate_of = primary.observation_id
    possible = defaultdict(list)
    for m in observations:
        if m.doi and m.inchikey and m.value_nm is not None:
            key = (m.doi.casefold(), m.inchikey, m.endpoint, m.relation, m.value_nm, m.cell_line.casefold())
            possible[key].append(m)
    for key, group in possible.items():
        if len({m.origin_assay_id for m in group}) > 1:
            for m in group:
                m.possible_duplicate_group = digest(key)[:16]
                m.warnings = sorted(set(m.warnings+['possible_cross_assay_duplicate_not_merged']))
    return observations


def export_run(config, store, rid):
    run = store.run(rid)
    observations = deduplicate(store.observations(rid))
    for m in observations:
        if m.extraction_method == "llm":
            if store.quality_passed(m.extractor_signature):
                m.review_reasons = [r for r in m.review_reasons if r != "extractor_not_validated"]
            elif "extractor_not_validated" not in m.review_reasons:
                m.review_reasons.append("extractor_not_validated")
    folder = config.root / "data/runs" / rid
    fields = list(Measurement.model_fields) + ["status"]
    counts = {}
    for name, rows in [("measurements", [m for m in observations if m.status == "accepted" and not m.duplicate_of]),
                       ("review", [m for m in observations if m.status == "review" and not m.duplicate_of]),
                       ("all_observations", observations)]:
        temp = folder / f"{name}.csv.tmp"
        with temp.open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for m in rows:
                row = m.model_dump()
                row["status"] = m.status
                for key,value in row.items():
                    if isinstance(value, (list,dict)):
                        row[key] = json.dumps(value, ensure_ascii=False)
                writer.writerow(row)
        temp.replace(folder / f"{name}.csv")
        counts[name] = len(rows)
    jobs = store.jobs(rid)
    report = {**run["report"], "run_id":rid, "status":run["status"], "virus":run["virus"],
              "inventory":run["inventory"], "options":run["options"], "exported_at":now(), "counts":counts,
              "unique_structures":len({m.inchikey for m in observations if m.inchikey}),
              "observed_extractor_signatures":dict(Counter(m.extractor_signature for m in observations)),
              "job_counts":dict(Counter(f"{j['source']}:{j['kind']}:{j['status']}" for j in jobs)),
              "review_reasons":dict(Counter(r for m in observations for r in m.review_reasons)),
              "unavailable_documents":[{"id":j["key"], "reason":"no_open_fulltext"} for j in jobs
                                       if j["kind"] == "paper" and (not j["payload"].get("pmcid") or j["payload"].get("isOpenAccess") != "Y")],
              "failed_jobs":[{"source":j["source"], "key":j["key"], "error":j["error"]} for j in jobs if j["status"] in ("failed","partial")],
              "quality_gate":{"literature_validated":store.quality_passed(config.extractor_signature),
                              "automatic_acceptance_is_not_manual_review":True},
              "limitations":["Coverage is restricted to the documented queries and connected sources.",
                             "Unprocessed papers and inaccessible supplements are not negative results.",
                             "No pooling across assays or automatic IC50/EC50/CC50 equivalence.",
                             "CC50 remains a separate measurement; no inferred activity/toxicity pairing."]}
    report["events"] = [dict(r) for r in store.db.execute("SELECT created,level,message FROM events WHERE run_id=?", (rid,))]
    (folder / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2))
    (folder / "manifest.json").write_text(json.dumps({"run":run, "config":config.settings,
                    "model":report.get("model", {}), "source_snapshots":sorted({m.source_snapshot for m in observations})},
                    ensure_ascii=False, indent=2))
    return folder, counts
