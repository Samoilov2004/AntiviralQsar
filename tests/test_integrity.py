from io import BytesIO
import zipfile
import pytest

from antiviral_qsar.chemistry import normalize, structure, validate_measurement
from antiviral_qsar.schema import Measurement
from antiviral_qsar.pipeline import deduplicate, known_context, export_run
from antiviral_qsar.documents import jats_chunks, attachment_chunks, zip_chunks
from antiviral_qsar.model import extract_json, grounded
from antiviral_qsar.sources import pubchem_unit


def measurement(**kwargs):
    defaults = dict(observation_id="1", source="test", source_record_id="1", source_url="https://example.org",
                    origin_assay_id="assay:1", smiles_original="CCO", endpoint="EC50", value_original="1.5",
                    units_original="uM", virus="Dengue virus", cell_line="Vero", evidence="a row",
                    evidence_grounded=True, cell_context_supported=True)
    return Measurement(**(defaults | kwargs))


@pytest.mark.parametrize("value,unit,relation,expected", [("1.5","uM","=",1500),("0.2","mM","=",200000),
    (">50","µM","=",50000),("1e-6","M","=",1000),("20","nM","<=",20),("250","pM","=",0.25)])
def test_unit_conversion_and_censoring(value, unit, relation, expected):
    nm, rel, p = normalize(value, unit, relation)
    assert nm == expected
    if value.startswith(">") or relation != "=":
        assert p is None
        assert rel != "="


@pytest.mark.parametrize("value,unit", [("NaN","uM"),("0","uM"),("-1","uM"),("1-3","uM"),
                                        ("1 ± 0.3","uM"),("5","ug/ml"),("70","%")])
def test_ambiguous_or_nonconcentration_values_are_not_labels(value, unit):
    with pytest.raises(ValueError):
        normalize(value, unit)


def test_relation_conflict_is_not_silently_resolved():
    with pytest.raises(ValueError, match="conflicting_relation"):
        normalize(">2", "uM", "<")


def test_stereochemistry_is_preserved_and_mixtures_rejected():
    a = structure("N[C@@H](C)C(=O)O")
    b = structure("N[C@H](C)C(=O)O")
    assert a[1] != b[1]
    assert structure("CC(=O)[O-].[Na+]")[2]
    with pytest.raises(ValueError, match="multiple_organic_fragments"):
        structure("CCO.CCOC")


def test_llm_results_need_both_evidence_and_quality_gate():
    m = validate_measurement(measurement(extraction_method="llm"))
    assert m.status == "review"
    assert "extractor_not_validated" in m.review_reasons
    m = validate_measurement(measurement(extraction_method="llm", evidence_grounded=False), llm_gate=True)
    assert "evidence_not_grounded" in m.review_reasons


def test_same_molecule_different_experiments_are_retained():
    a = validate_measurement(measurement())
    b = validate_measurement(measurement(observation_id="2", origin_assay_id="assay:2"))
    result = deduplicate([a,b])
    assert not any(m.duplicate_of for m in result)


def test_imported_copy_merges_provenance_not_independent_measurements():
    a = validate_measurement(measurement(source="chembl"))
    b = validate_measurement(measurement(observation_id="2", source="pubchem"))
    result = deduplicate([a,b])
    assert sum(bool(m.duplicate_of) for m in result) == 1
    assert len(a.provenance) == 2


def test_organism_label_does_not_make_biochemical_assay_cellular():
    result = known_context("Dissociation rate constant for Dengue virus envelope protein binding by surface plasmon resonance")
    assert result.kind == "biochemical"
    result = known_context("Antiviral activity against Dengue virus in infected Vero cells")
    assert result.kind == "cellular_antiviral" and result.cell_line == "Vero"


def test_table_units_headers_and_footnotes_survive_chunking():
    xml = b'''<article><front><article-title>Antiviral compounds</article-title></front><body>
    <p>Dengue virus was tested in Vero cells.</p><table-wrap id="T1"><caption><p>EC50 in uM</p></caption>
    <table><tr><th rowspan="2">Compound</th><th colspan="2">Vero</th></tr>
    <tr><th>EC50</th><th>CC50</th></tr><tr><td>7b</td><td>1.2</td><td>&gt;50</td></tr></table>
    <table-wrap-foot><p>Exposure 48 hours.</p></table-wrap-foot></table-wrap></body></article>'''
    chunks, _ = jats_chunks(xml)
    table = chunks[0]
    assert "EC50 in uM" in table.text and "Exposure 48 hours" in table.text
    assert "7b | 1.2 | >50" in table.text
    assert "Vero cells" in table.context


def test_empty_scan_and_zip_path_are_reported_not_executed(tmp_path):
    import pymupdf
    document = pymupdf.open()
    document.new_page()
    chunks, warnings = attachment_chunks("scan.pdf", document.tobytes())
    assert not chunks and "image_only_page" in warnings[0]
    output = BytesIO()
    with zipfile.ZipFile(output, "w") as z:
        z.writestr("../../evil.sh", "touch bad")
        z.writestr("values.csv", "Compound,EC50\n7b,1.2")
    parsed, issues = zip_chunks(output.getvalue())
    assert parsed and any("unsupported_attachment" in w for w in issues)
    assert not (tmp_path / "evil.sh").exists()


def test_model_evidence_is_checked_and_incomplete_json_rejected():
    assert grounded("IC50 1.2 µM", "Measured IC50   1.2 μM today")
    assert not grounded("IC50 5", "IC50 1.2")
    assert extract_json('```json\n{"measurements":[]}\n```') == {"measurements":[]}
    with pytest.raises(ValueError):
        extract_json('{"measurements": [')


def test_pubchem_unit_enum_preserves_nanomolar_vs_micromolar():
    assert pubchem_unit(5) == "uM"
    assert pubchem_unit(6) == "nM"


def test_resume_retains_done_jobs_and_does_not_duplicate_rows(store, config):
    rid = store.new_run({"name":"Dengue virus"}, {})
    store.add_job(rid,"chembl","assay","A",{})
    store.add_job(rid,"chembl","assay","B",{})
    jobs = store.jobs(rid)
    store.start_job(jobs[0])
    store.finish_job(jobs[0],[validate_measurement(measurement())])
    store.start_job(jobs[1])
    assert [j["key"] for j in store.jobs(rid,unfinished=True)] == ["B"]
    store.finish_job(jobs[0],[validate_measurement(measurement())])
    assert len(store.observations(rid)) == 1
    folder,counts = export_run(config,store,rid)
    assert counts["measurements"] == 1
    assert (folder/"report.json").exists()
