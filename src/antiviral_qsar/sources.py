from __future__ import annotations

import re
from urllib.parse import quote
from lxml import etree

from .network import Network
from .schema import Virus

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"


def clean_text(element):
    return " ".join(" ".join(element.itertext()).split())


def xml_parse(data):
    return etree.fromstring(data, parser=etree.XMLParser(resolve_entities=False, no_network=True, recover=False))


def search_terms(virus: Virus):
    result = []
    for term in [virus.query, virus.name, *virus.synonyms]:
        term = term.strip()
        if len(term) >= 3 and term.lower() not in {x.lower() for x in result}:
            result.append(term)
    return result


def taxonomy_candidates(net: Network, query: str, tax_id: str | None = None):
    if tax_id:
        ids = [tax_id]
    else:
        result, _ = net.json(f"{EUTILS}/esearch.fcgi", {"db":"taxonomy", "term":f'"{query}"[All Names]',
                                                       "retmode":"json", "retmax":30}, max_age=86400)
        ids = result.get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise ValueError("Name not found in NCBI Taxonomy. Use an English name or --tax-id.")
    data, _ = net.get(f"{EUTILS}/efetch.fcgi", {"db":"taxonomy", "id":",".join(ids), "retmode":"xml"})
    root = xml_parse(data)
    candidates = []
    for taxon in root.findall("Taxon"):
        ancestors = [x.text for x in taxon.findall("LineageEx/Taxon/TaxId")]
        if "10239" not in ancestors and taxon.findtext("TaxId") != "10239":
            continue
        synonyms = []
        for node in taxon.findall("OtherNames/*"):
            if node.tag in ("Synonym", "EquivalentName", "GenbankCommonName", "CommonName"):
                name = clean_text(node)
                if name:
                    synonyms.append(name)
        candidates.append(Virus(query=query, name=taxon.findtext("ScientificName", ""),
                                tax_id=taxon.findtext("TaxId", ""), synonyms=sorted(set(synonyms))))
    if not candidates:
        raise ValueError("No viral taxon with this name was found in NCBI Taxonomy")
    exact = [v for v in candidates if query.casefold() == v.name.casefold()]
    return exact if len(exact) == 1 else candidates


class ChEMBL:
    name = "chembl"

    def __init__(self, net):
        self.net = net
        self._molecules = {}

    def pages(self, resource, params, *, fresh=False):
        url = f"{CHEMBL}/{resource}.json"
        query = {**params, "limit":1000}
        while url:
            data, snapshot = self.net.json(url, query, max_age=86400 if fresh else None)
            key = {"activity":"activities", "assay":"assays", "molecule":"molecules",
                   "compound_record":"compound_records", "document":"documents"}.get(resource, resource+"s")
            for row in data.get(key, []):
                yield row, snapshot
            next_url = data.get("page_meta", {}).get("next")
            url = ("https://www.ebi.ac.uk" + next_url if next_url and next_url.startswith("/") else next_url)
            query = None

    def discover(self, virus):
        status, _ = self.net.json(f"{CHEMBL}/status.json", max_age=86400)
        seen, assays = set(), []
        queries = [{"assay_tax_id":virus.tax_id}]
        queries += [{"description__icontains":term} for term in search_terms(virus)]
        for query in queries:
            for row, snapshot in self.pages("assay", query, fresh=True):
                key = row["assay_chembl_id"]
                if key not in seen:
                    seen.add(key)
                    row["_snapshot"] = snapshot
                    assays.append(row)
        return assays, {"status":"available", "release":status.get("chembl_db_version"),
                        "query_count":len(queries), "assays":len(assays)}

    def activities(self, assay_id):
        return self.pages("activity", {"assay_chembl_id":assay_id, "standard_type__in":"IC50,EC50,CC50"})

    def molecule(self, compound_id):
        if compound_id in self._molecules:
            return self._molecules[compound_id]
        data, _ = self.net.json(f"{CHEMBL}/molecule/{compound_id}.json")
        self._molecules[compound_id] = data
        return data

    def prefetch_molecules(self, ids):
        ids = sorted(set(ids) - self._molecules.keys())
        for start in range(0, len(ids), 100):
            for row, _ in self.pages("molecule", {"molecule_chembl_id__in":','.join(ids[start:start+100])}):
                self._molecules[row['molecule_chembl_id']] = row

    def document(self, doc_id):
        if not doc_id:
            return {}
        data, _ = self.net.json(f"{CHEMBL}/document/{doc_id}.json")
        return data

    def documents_by_reference(self, doi="", pmid=""):
        query = {"doi__iexact":doi} if doi else {"pubmed_id":pmid} if pmid else None
        return [r for r, _ in self.pages("document", query)] if query else []

    def linked_cytotoxicity(self, doc_ids):
        seen = set()
        for doc_id in sorted(set(doc_ids)):
            for assay, snapshot in self.pages("assay", {"document_chembl_id":doc_id}):
                description = assay.get("description", "").lower()
                if re.search(r"cytotoxi|cell viability|cc50|tc50", description):
                    aid = assay["assay_chembl_id"]
                    if aid not in seen:
                        seen.add(aid)
                        assay["_snapshot"] = snapshot
                        yield assay

    def compound_in_document(self, doc_id, label):
        if not doc_id:
            return []
        found = []
        for record, _ in self.pages("compound_record", {"document_chembl_id":doc_id}):
            if str(record.get("compound_name", "")).strip().casefold() == label.strip().casefold():
                found.append(record["molecule_chembl_id"])
        return sorted(set(found))


class PubChem:
    name = "pubchem"

    def __init__(self, net):
        self.net = net
        self._substances = {}
        self._compounds = {}

    def discover(self, virus):
        quoted = [f'"{t.replace(chr(34), "")}"[All Fields]' for t in search_terms(virus)]
        term = " OR ".join(quoted + [f"{virus.tax_id}[Taxonomy ID]"])
        result, _ = self.net.json(f"{EUTILS}/esearch.fcgi", {"db":"pcassay", "term":term,
                                   "retmode":"json", "retmax":10000}, max_age=86400)
        response = result["esearchresult"]
        if response.get("ERROR"):
            raise ValueError(str(response["ERROR"]))
        ids = response.get("idlist", [])
        total = int(response.get("count", 0))
        for offset in range(10000, total, 10000):
            more, _ = self.net.json(f"{EUTILS}/esearch.fcgi", {"db":"pcassay", "term":term,
                      "retmode":"json", "retmax":10000, "retstart":offset}, max_age=86400)
            ids.extend(more["esearchresult"]["idlist"])
        ids = list(dict.fromkeys(ids))
        if len(ids) != total:
            raise ValueError(f"Incomplete PubChem inventory: {len(ids)} / {total}")
        return ids, {"status":"available", "assays":len(ids), "query":term}

    def description(self, aid):
        data, snapshot = self.net.json(f"{PUBCHEM}/assay/aid/{aid}/description/JSON")
        return self.containers(data)[0]["assay"]["descr"], snapshot

    @staticmethod
    def containers(data):
        if "PC_AssaySubmit" in data:
            return [data["PC_AssaySubmit"]]
        if "PC_AssayContainer" in data:
            return data["PC_AssayContainer"]
        raise ValueError("unrecognized_pubchem_assay_response")

    def results(self, aid):
        # Request by SID batches: full-assay PUG REST is limited to 10,000 rows.
        info, _ = self.net.json(f"{PUBCHEM}/assay/aid/{aid}/sids/JSON")
        entries = info.get("InformationList", {}).get("Information", [])
        sids = [str(sid) for entry in entries for sid in entry.get("SID", [])]
        if not sids:
            return
        for start in range(0, len(sids), 200):
            batch = sids[start:start+200]
            data, snapshot = self.net.json(f"{PUBCHEM}/assay/aid/{aid}/JSON", {"sid":",".join(batch)})
            seen = set()
            for container in self.containers(data):
                for row in container.get("data", []):
                    seen.add(str(row["sid"]))
                    yield row, snapshot
            if seen != set(batch):
                raise ValueError(f"incomplete_pubchem_assay_results:{len(seen)}/{len(batch)}")

    def sid_compound(self, sid):
        if str(sid) in self._substances:
            ids = self._substances[str(sid)]
            return (self.compound(ids[0]), "") if len(ids) == 1 else ({}, "ambiguous_substance_structure")
        info, _ = self.net.json(f"{PUBCHEM}/substance/sid/{sid}/cids/JSON", {"cids_type":"standardized"})
        entries = info.get("InformationList", {}).get("Information", [])
        ids = sorted({cid for e in entries for cid in e.get("CID", [])})
        if len(ids) != 1:
            return {}, "ambiguous_substance_structure"
        return self.compound(ids[0]), ""

    def compound(self, cid):
        if int(cid) in self._compounds:
            return self._compounds[int(cid)]
        data, _ = self.net.json(f"{PUBCHEM}/compound/cid/{cid}/property/SMILES,InChIKey/JSON")
        row = data["PropertyTable"]["Properties"][0]
        self._compounds[int(cid)] = row
        return row

    def prefetch_substances(self, sids):
        sids = sorted(set(map(str, sids)) - self._substances.keys())
        for start in range(0, len(sids), 100):
            batch = ','.join(sids[start:start+100])
            data, _ = self.net.json(f"{PUBCHEM}/substance/sid/{batch}/cids/JSON", {"cids_type":"standardized"})
            for row in data.get('InformationList', {}).get('Information', []):
                self._substances[str(row['SID'])] = sorted(set(row.get('CID', [])))
        ids = sorted({cid for values in self._substances.values() if len(values) == 1 for cid in values}
                     - self._compounds.keys())
        for start in range(0, len(ids), 100):
            batch = ','.join(map(str, ids[start:start+100]))
            data, _ = self.net.json(f"{PUBCHEM}/compound/cid/{batch}/property/SMILES,InChIKey/JSON")
            for row in data['PropertyTable']['Properties']:
                self._compounds[int(row['CID'])] = row

    def resolve_name(self, name):
        if re.fullmatch(r"(?i)(?:compound\s*)?\d+[a-z]*", name.strip()):
            return {}, "document_local_compound_label"
        info, _ = self.net.json(f"{PUBCHEM}/compound/name/{quote(name, safe='')}/cids/JSON", {"name_type":"complete"})
        ids = info.get("IdentifierList", {}).get("CID", [])
        if len(ids) != 1:
            return {}, "ambiguous_compound_name"
        cid = ids[0]
        synonyms, _ = self.net.json(f"{PUBCHEM}/compound/cid/{cid}/synonyms/JSON")
        names = synonyms.get("InformationList", {}).get("Information", [{}])[0].get("Synonym", [])
        if name.strip().casefold() not in {n.strip().casefold() for n in names}:
            return {}, "compound_name_not_exact_synonym"
        return self.compound(cid), ""


class EuropePMC:
    name = "europepmc"

    def __init__(self, net):
        self.net = net

    def discover(self, virus, references=()):
        terms = " OR ".join('"'+t.replace('"', '')+'"' for t in search_terms(virus))
        query = f"({terms}) AND (IC50 OR EC50 OR CC50 OR antiviral OR cytotoxicity)"
        cursor, seen, results = "*", set(), []
        queries = [query]
        for start in range(0, len(references), 50):
            batch = references[start:start+50]
            if batch:
                queries.append(" OR ".join(f'EXT_ID:{r}' for r in batch))
        hit_count = 0
        for search in queries:
            cursor = "*"
            while cursor:
                data, _ = self.net.json(f"{EPMC}/search", {"query":search, "format":"json",
                      "pageSize":1000, "cursorMark":cursor, "resultType":"core"}, max_age=86400)
                if search == query:
                    hit_count = data.get("hitCount", 0)
                batch = data.get("resultList", {}).get("result", [])
                for row in batch:
                    key = paper_identity(row)
                    if key not in seen:
                        seen.add(key)
                        results.append(row)
                next_cursor = data.get("nextCursorMark")
                if not batch or next_cursor == cursor:
                    break
                cursor = next_cursor
        return results, {"status":"available", "query":query, "search_hits":hit_count,
                         "publications":len(results)}

    def fulltext(self, paper):
        pmcid = paper.get("pmcid")
        if not pmcid or paper.get("isOpenAccess") != "Y":
            raise ValueError("open_fulltext_unavailable")
        try:
            return self.net.get(f"{EPMC}/{pmcid}/fullTextXML")
        except Exception as xml_error:
            # Use only publisher/PMC PDF links declared by the open-access metadata.
            for link in paper.get('fullTextUrlList', {}).get('fullTextUrl', []):
                if link.get('documentStyle', '').lower() != 'pdf' or link.get('availability', '').lower() != 'open access':
                    continue
                try:
                    raw, snapshot = self.net.get(link['url'], max_bytes=50_000_000)
                    if raw.lstrip().startswith(b'%PDF-'):
                        return raw, snapshot
                except Exception:
                    continue
            raise ValueError(f"open_fulltext_unavailable:{xml_error}") from xml_error

    def supplements(self, paper):
        pmcid = paper.get("pmcid")
        if not pmcid:
            return None, ""
        return self.net.get(f"{EPMC}/{pmcid}/supplementaryFiles", max_bytes=100_000_000)


def paper_identity(paper):
    doi = paper.get("doi", "").strip().lower().removeprefix("https://doi.org/")
    return "doi:" + doi if doi else paper.get("pmcid") or f"{paper.get('source')}:{paper.get('id')}"


def pubchem_value(field):
    value = field.get("value", {})
    for kind in ("fval", "ival", "sval", "bval"):
        if kind in value:
            return str(value[kind])
    return ""


def pubchem_unit(code):
    # NCBI official pcassay.asn, PC-ResultType.unit (not the TID number).
    return {4:"mM", 5:"uM", 6:"nM", 7:"pM", 8:"fM", 9:"mg/mL", 10:"ug/mL",
            11:"ng/mL", 12:"pg/mL", 13:"fg/mL", 14:"M", 15:"%"}.get(code, str(code or ""))
