"""Publication coverage in the current database collection, before local inference."""
from __future__ import annotations

from collections import defaultdict
import re
from urllib.parse import unquote


def normalize_doi(value):
    value = unquote(str(value or "").strip())
    value = re.sub(r"^(?:doi:\s*|https?://(?:dx\.)?doi\.org/)", "", value, flags=re.I)
    value = value.strip().casefold()
    return value if re.fullmatch(r"10\.\d{4,9}/\S+", value) else ""


def normalize_pmid(value):
    value = re.sub(r"^pmid:\s*", "", str(value or "").strip(), flags=re.I)
    return str(int(value)) if value.isascii() and value.isdigit() and int(value) > 0 else ""


def reference_keys(doi="", pmids=()):
    keys = {"pmid:" + p for value in pmids if (p := normalize_pmid(value))}
    if value := normalize_doi(doi):
        keys.add("doi:" + value)
    return keys


def paper_pmids(paper):
    # Europe PMC IDs from PMC, patents, etc. are not PubMed IDs.
    value = paper.get("pmid") or (paper.get("id") if paper.get("source") == "MED" else "")
    return [value] if normalize_pmid(value) else []


def measurement_keys(row):
    legacy_pmids = re.findall(r"(?:^|;)PMID:(\d+)(?=;|$)", row.document_id, re.I)
    return reference_keys(row.doi, [*row.pmids, *legacy_pmids])


def conflicting_reference(keys, reference):
    if not keys.intersection(reference):
        return False
    for prefix in ('doi:', 'pmid:'):
        requested = {key for key in keys if key.startswith(prefix)}
        stored = {key for key in reference if key.startswith(prefix)}
        if requested and stored and requested != stored:
            return True
    return False


def missing_fields(rows):
    if not rows:
        return ["no_imported_measurements"]
    gaps, endpoints = set(), defaultdict(set)
    for row in rows:
        for field in ("smiles", "cell_line", "duration", "method"):
            if not getattr(row, field).strip():
                gaps.add("missing_" + field)
        if row.endpoint != "CC50":
            for field in ("virus", "strain"):
                if not getattr(row, field).strip():
                    gaps.add("missing_" + field)
        if row.value_nm is None:
            gaps.add("missing_usable_value")
        if not row.cell_context_supported:
            gaps.add("unverified_cellular_context")
        compound = row.inchikey or row.compound_id or row.observation_id
        endpoints[compound].add(row.endpoint)
    for present in endpoints.values():
        if "CC50" not in present:
            gaps.add("missing_CC50")
        if not present.intersection({"IC50", "EC50"}):
            gaps.add("missing_antiviral_endpoint")
    return sorted(gaps)


class PublicationCoverage:
    def __init__(self, jobs, observations, *, incomplete=False):
        self.jobs = {}
        self.rows = {m.observation_id: m for m in observations if m.source in ("chembl", "pubchem")}
        self.assays_by_ref, self.rows_by_ref = defaultdict(set), defaultdict(set)
        self.ambiguous_rows = set()
        self.incomplete = incomplete
        by_assay, by_document = {}, defaultdict(list)
        for job in jobs:
            if job["source"] not in ("chembl", "pubchem") or job["kind"] not in ("assay", "cytotoxicity"):
                continue
            key = job["source"] + ":" + job["key"]
            self.jobs[key] = job
            by_assay[(job["source"], job["key"])] = job
            if doc := job["payload"].get("document_chembl_id"):
                by_document[doc].append(job)
            refs = job["payload"].get("publication_references", [])
            for ref in refs:
                for token in reference_keys(ref.get("doi"), [ref.get("pmid")]):
                    self.assays_by_ref[token].add(key)
            if job["payload"].get("publication_reference_status") != "available":
                self.incomplete = True
        for row in self.rows.values():
            keys = measurement_keys(row)
            source_assay = row.source_assay_id or row.source_url.rstrip("/").rsplit("/", 1)[-1]
            linked = [by_assay[(row.source, source_assay)]] if (row.source, source_assay) in by_assay else []
            if row.source == "chembl":
                linked += by_document.get(row.document_id, [])
            for job in linked:
                refs = job["payload"].get("publication_references", [])
                for ref in refs:
                    keys |= reference_keys(ref.get("doi"), [ref.get("pmid")])
                if len(refs) > 1:
                    self.ambiguous_rows.add(row.observation_id)
            if sum(k.startswith("pmid:") for k in keys) > 1:
                self.ambiguous_rows.add(row.observation_id)
            for token in keys:
                self.rows_by_ref[token].add(row.observation_id)

    def plan(self, paper):
        keys = reference_keys(paper.get("doi"), paper_pmids(paper))
        matched = {key for key in keys if key in self.assays_by_ref or key in self.rows_by_ref}
        assays = set().union(*(self.assays_by_ref[k] for k in matched))
        ids = set().union(*(self.rows_by_ref[k] for k in matched))
        gaps = []
        if matched:
            gaps = missing_fields([self.rows[i] for i in sorted(ids)])
            references = [reference_keys(ref.get('doi'), [ref.get('pmid')])
                          for a in assays for ref in self.jobs[a]['payload'].get('publication_references', [])]
            references += [measurement_keys(self.rows[i]) for i in ids if i not in self.ambiguous_rows]
            if any(conflicting_reference(keys, ref) for ref in references):
                gaps.append('conflicting_publication_identifiers')
            if ids & self.ambiguous_rows or any(
                len(self.jobs[a]["payload"].get("publication_references", [])) > 1 for a in assays
            ):
                gaps.append("ambiguous_database_references")
            if any(self.jobs[a]["status"] in ("failed", "partial", "pending", "running") for a in assays):
                gaps.append("database_assay_incomplete")
            if self.incomplete:
                gaps.append("database_reference_coverage_incomplete")
            mode = "supplement" if gaps else "skip"
        else:
            mode = "new" if keys and not self.incomplete else "unmatched"
        return {"mode": mode, "identifiers": sorted(keys), "matched_identifiers": sorted(matched),
                "database_assays": sorted(assays), "database_observations": sorted(ids),
                "gaps": sorted(set(gaps)), "reference_coverage_incomplete": self.incomplete,
                "reason": "no_detected_gaps" if mode == "skip" else
                          "database_gaps" if mode == "supplement" else
                          "no_matching_database_reference" if mode == "new" else
                          "missing_identifiers_or_incomplete_reference_coverage"}

    def supplement(self, records, plan):
        """Retain raw extractions; avoid re-exporting information already imported.

        A covered row is not asserted to be the same experiment. It remains in
        all_observations.csv. Missing fields are never copied into database rows.
        """
        conflict = 'conflicting_publication_identifiers' in plan.get('gaps', [])
        baseline = [self.rows[i] for i in plan["database_observations"] if i not in self.ambiguous_rows and not conflict]
        for row in records:
            row.literature_role, row.supplement_of = "", []
            row.review_reasons = [r for r in row.review_reasons if r not in (
                "supplement_link_requires_review", "publication_identity_requires_review")]
            if conflict:
                row.review_reasons.append('publication_identity_requires_review')
            if plan["mode"] in ("new", "unmatched"):
                continue
            row.literature_role = "additional_measurement"
            candidates = [m for m in baseline if same_compound(row, m) and row.endpoint == m.endpoint]
            compatible = [m for m in candidates if compatible_context(row, m)]
            same_value = [m for m in compatible if row.value_nm is not None and
                          row.value_nm == m.value_nm and row.relation == m.relation]
            if any(covers(m, row) for m in same_value):
                row.literature_role = "already_covered"
                row.supplement_of = [m.observation_id for m in same_value if covers(m, row)]
            elif same_value and any(any(getattr(row, f) and not getattr(m, f)
                                       for f in ('smiles', *context_fields(row))) for m in same_value):
                row.literature_role = "field_supplement"
                row.supplement_of = [m.observation_id for m in same_value]
                row.review_reasons.append("supplement_link_requires_review")
            elif compatible:
                # A different number or inequality may be another experiment.
                # Preserve it for review; never replace the structured result.
                row.literature_role = "unresolved_overlap"
                row.supplement_of = [m.observation_id for m in compatible]
                row.review_reasons.append("supplement_link_requires_review")
        return records


def same_compound(a, b):
    if a.inchikey and b.inchikey:
        return a.inchikey == b.inchikey
    return bool(a.compound_id and a.compound_id == b.compound_id)


def text_key(value):
    return " ".join(value.split()).casefold()


def context_fields(row):
    return ("cell_line", "duration", "method") + (() if row.endpoint == "CC50" else ("virus", "strain"))


def compatible_context(a, b):
    return all(not getattr(a, f) or not getattr(b, f) or text_key(getattr(a, f)) == text_key(getattr(b, f))
               for f in context_fields(a))


def covers(database, extracted):
    # Do not hide an extraction with unresolved structure/context, or use an
    # uncertain database row as evidence that a measurement is covered.
    if database.review_reasons or not extracted.inchikey or not extracted.cell_line:
        return False
    return all(not getattr(extracted, f) or text_key(getattr(database, f)) == text_key(getattr(extracted, f))
               for f in context_fields(extracted))
