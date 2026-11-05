import csv
import json

import pytest

from antiviral_qsar.chemistry import validate_measurement
from antiviral_qsar.literature import PublicationCoverage, normalize_doi, normalize_pmid, paper_pmids
from antiviral_qsar.pipeline import Pipeline, export_run
from antiviral_qsar.schema import Virus
from test_integrity import measurement


DOI = '10.1234/vacv.1'
VIRUS = Virus(query='Vaccinia virus', name='Vaccinia virus', tax_id='10245')


def record(identifier='ec', **changes):
    return validate_measurement(measurement(**(dict(
        observation_id=identifier, source='chembl', source_record_id=identifier, source_assay_id='A',
        source_url='https://www.ebi.ac.uk/chembl/explore/assay/A', document_id='D', doi=DOI, pmids=['123'],
        virus=VIRUS.name, strain='WR', cell_line='Vero', duration='48 h', method='plaque reduction',
        compound_id='CHEMBL1', endpoint='EC50') | changes)))


def complete_rows():
    return [record(), record('cc', endpoint='CC50', value_original='50', virus='', strain='')]


def paper(identifier='123', **changes):
    return dict(source='MED', id=identifier, pmcid='PMC' + identifier, doi=DOI,
                isOpenAccess='Y', title='Vaccinia virus IC50 EC50 CC50', **changes)


def database_job(refs=None, status='done', source='chembl'):
    return dict(source=source, kind='assay', key='A', status=status, payload={
        'document_chembl_id': 'D' if source == 'chembl' else '',
        'publication_references': refs if refs is not None else [{'doi': DOI, 'pmid': '123'}],
        'publication_reference_status': 'available'})


def extracted(identifier='lit', **changes):
    return record(identifier, **(dict(source='europepmc', source_assay_id='', document_id='PMC123',
                                     origin_assay_id='paper:table-1', extraction_method='llm') | changes))


def test_reference_normalization_doi_or_pmid_and_no_pmc_id_as_pmid():
    assert normalize_doi(' HTTPS://DX.DOI.ORG/10.1234/VACV%2E1 ') == DOI
    assert normalize_pmid('PMID: 00123') == '123'
    assert normalize_pmid('PMC123') == ''
    assert paper_pmids({'source': 'PMC', 'id': '123'}) == []
    assert normalize_doi('') == ''
    coverage = PublicationCoverage([database_job()], complete_rows())
    for candidate in ({'doi': 'DOI:10.1234/VACV.1'}, {'source': 'MED', 'id': '123'}):
        assert coverage.plan(candidate)['mode'] == 'skip'
    assert coverage.plan({'source': 'PMC', 'id': '123'})['mode'] == 'unmatched'
    assert coverage.plan({'doi': '10.1234/other'})['mode'] == 'new'


def test_old_pubchem_pmid_list_matches_without_new_schema_fields():
    row = record(source='pubchem', doi='', pmids=[], document_id='PMID:123', source_assay_id='')
    plan = PublicationCoverage([], [row]).plan({'source': 'MED', 'id': '123'})
    assert plan['mode'] == 'supplement'
    assert plan['database_observations'] == ['ec']


def test_old_chembl_rows_link_by_document_when_doi_and_pmid_missing():
    rows = [m.model_copy(update={'doi': '', 'pmids': [], 'source_assay_id': '',
                                'source_url': 'https://example.org/old'}) for m in complete_rows()]
    assert PublicationCoverage([database_job()], rows).plan({'pmid': '123'})['mode'] == 'skip'


@pytest.mark.parametrize('missing,reason', [('cell_line', 'missing_cell_line'), ('smiles', 'missing_smiles'),
                                         ('duration', 'missing_duration'), ('strain', 'missing_strain')])
def test_known_articles_with_field_gaps_are_eligible(missing, reason):
    rows = complete_rows()
    setattr(rows[0], missing, '')
    plan = PublicationCoverage([database_job()], rows).plan(paper())
    assert plan['mode'] == 'supplement' and reason in plan['gaps']


def test_missing_cc50_is_a_gap_but_missing_ic50_with_ec50_is_not():
    plan = PublicationCoverage([database_job()], [record()]).plan(paper())
    assert plan['gaps'] == ['missing_CC50']
    plan = PublicationCoverage([database_job()], complete_rows()).plan(paper())
    assert plan['gaps'] == []
    plan = PublicationCoverage([database_job()], []).plan(paper())
    assert plan['gaps'] == ['no_imported_measurements']


def test_failed_assays_and_unknown_reference_coverage_never_skip_known_papers():
    plan = PublicationCoverage([database_job(status='failed')], complete_rows()).plan(paper())
    assert 'database_assay_incomplete' in plan['gaps']
    plan = PublicationCoverage([database_job()], complete_rows(), incomplete=True).plan(paper())
    assert plan['mode'] == 'supplement'
    assert PublicationCoverage([], [], incomplete=True).plan({'doi': DOI})['mode'] == 'unmatched'


def test_multiple_pubmed_references_do_not_imply_complete_coverage_of_each_article():
    job = database_job([{'pmid': '123'}, {'pmid': '456'}], source='pubchem')
    rows = [m.model_copy(update={'source': 'pubchem', 'doi': '', 'pmids': ['123', '456']})
            for m in complete_rows()]
    coverage = PublicationCoverage([job], rows)
    plan = coverage.plan({'pmid': '123'})
    assert plan['mode'] == 'supplement'
    assert 'ambiguous_database_references' in plan['gaps']
    assert coverage.supplement([extracted()], plan)[0].literature_role == 'additional_measurement'


def test_conflicting_doi_and_pmid_cannot_hide_a_paper_or_its_extractions():
    coverage = PublicationCoverage([database_job()], complete_rows())
    plan = coverage.plan({'doi': DOI, 'pmid': '456'})
    assert plan['mode'] == 'supplement' and 'conflicting_publication_identifiers' in plan['gaps']
    row = coverage.supplement([extracted()], plan)[0]
    assert row.literature_role != 'already_covered'
    assert 'publication_identity_requires_review' in row.review_reasons


def test_supplement_retains_changes_and_conflicts_without_overwriting_database():
    rows = complete_rows()
    rows[0].duration = ''
    coverage = PublicationCoverage([database_job()], rows)
    plan = coverage.plan(paper())
    candidates = [extracted(), extracted('covered', endpoint='CC50', value_original='50', virus='', strain=''),
                  extracted('new-cell', cell_line='HeLa'), extracted('other-value', value_original='2'),
                  extracted('other-sign', relation='>'), extracted('other-endpoint', endpoint='IC50')]
    result = coverage.supplement(candidates, plan)
    assert [m.literature_role for m in result] == ['field_supplement', 'already_covered',
        'additional_measurement', 'unresolved_overlap', 'unresolved_overlap', 'additional_measurement']
    assert result[0].supplement_of == ['ec']
    assert 'supplement_link_requires_review' in result[0].review_reasons
    assert rows[0].duration == '' and rows[0].value_original == '1.5'
    assert result[4].relation == '>'


def test_covered_rows_survive_audit_export_and_cached_rows_are_reassessed(config, store):
    coverage = PublicationCoverage([database_job()], complete_rows())
    lit = coverage.supplement([extracted()], coverage.plan(paper()))[0]
    assert lit.literature_role == 'already_covered'
    rid = store.new_run(VIRUS.model_dump(), {})
    store.add_job(rid, 'europepmc', 'paper', 'p', paper())
    store.finish_job(store.jobs(rid)[0], [*complete_rows(), lit])
    folder, counts = export_run(config, store, rid)
    assert counts == {'measurements': 2, 'review': 0, 'all_observations': 3}
    with (folder / 'all_observations.csv').open(encoding='utf-8-sig') as f:
        assert any(r['literature_role'] == 'already_covered' for r in csv.DictReader(f))
    revised = coverage.supplement([lit], {'mode': 'new', 'database_observations': []})[0]
    assert not revised.literature_role and revised.supplement_of == []


def setup_pipeline(config, store, monkeypatch, rows, papers):
    rid = store.new_run(VIRUS.model_dump(), {'papers': 'all'})
    store.add_job(rid, 'chembl', 'assay', 'A', {'document_chembl_id': 'D'})
    store.finish_job(store.jobs(rid)[0], rows)
    store.add_jobs(rid, 'europepmc', 'paper', papers)
    pipeline = Pipeline(config, store, use_model=False)
    monkeypatch.setattr(pipeline.chembl, 'document', lambda did: {'doi': DOI, 'pubmed_id': 123})
    monkeypatch.setattr(pipeline.chembl, 'linked_cytotoxicity', lambda ids: [])
    return rid, pipeline


def test_unmatched_papers_take_budget_before_highly_ranked_database_papers(config, store, monkeypatch):
    known = paper()
    new = dict(paper('456'), doi='10.1234/new', title='Short title')
    rid, pipeline = setup_pipeline(config, store, monkeypatch, [record()], [('known', known), ('new', new)])
    calls = []
    monkeypatch.setattr(pipeline, 'paper', lambda p, v, r: (calls.append(p['pmcid']) or [], []))
    pipeline.execute(rid, '1')
    assert calls == ['PMC456']
    assert store.run(rid)['options']['selected_paper_keys'] == ['new']
    pipeline.execute(rid, '1')
    assert calls == ['PMC456', 'PMC123']
    pipeline.net.close()


def test_complete_database_paper_is_skipped_before_downloading_and_budget_counts_only_selected(config, store, monkeypatch):
    new = dict(paper('456'), doi='10.1234/new')
    rid, pipeline = setup_pipeline(config, store, monkeypatch, complete_rows(), [('known', paper()), ('new', new)])
    calls = []
    monkeypatch.setattr(pipeline, 'paper', lambda p, v, r: (calls.append(p['pmcid']) or [], []))
    folder, _ = pipeline.execute(rid, '1')
    assert calls == ['PMC456']
    known = next(j for j in store.jobs(rid, 'paper') if j['key'] == 'known')
    assert known['status'] == 'skipped' and known['attempts'] == 0
    assert known['payload']['collection_plan']['reason'] == 'no_detected_gaps'
    report = json.loads((folder / 'report.json').read_text())
    assert report['paper_coverage']['counts'] == {'skip': 1, 'new': 1}
    pipeline.execute(rid, None)
    assert calls == ['PMC456']
    pipeline.net.close()


def test_reference_lookup_failure_is_reported_and_does_not_suppress_paper(config, store, monkeypatch):
    rid, pipeline = setup_pipeline(config, store, monkeypatch, complete_rows(), [('known', paper())])
    def fail(did):
        raise RuntimeError('source unavailable')
    monkeypatch.setattr(pipeline.chembl, 'document', fail)
    calls = []
    monkeypatch.setattr(pipeline, 'paper', lambda p, v, r: (calls.append(p['pmcid']) or [], []))
    folder, _ = pipeline.execute(rid, 'all')
    assert calls == ['PMC123']
    report = json.loads((folder / 'report.json').read_text())
    assert report['paper_coverage']['reference_errors'][0]['assay'] == 'A'
    assert report['paper_coverage']['papers'][0]['mode'] == 'supplement'
    pipeline.net.close()


def test_pubchem_reference_lookup_handles_nested_and_direct_pmids(config, store, monkeypatch):
    rid = store.new_run(VIRUS.model_dump(), {})
    store.add_job(rid, 'pubchem', 'assay', '42', {})
    store.finish_job(store.jobs(rid)[0])
    pipeline = Pipeline(config, store, use_model=False)
    monkeypatch.setattr(pipeline.pubchem, 'description', lambda aid: (
        {'xref': [{'xref': {'pmid': 123}}, {'pmid': 456}, {'pmid': 123}, {'taxonomy': 10245}]}, 'raw'))
    coverage = pipeline.publication_coverage(rid)
    assert coverage.plan({'pmid': '456'})['database_assays'] == ['pubchem:42']
    assert len(store.jobs(rid)[0]['payload']['publication_references']) == 2
    pipeline.net.close()


def test_cached_paper_is_reconciled_even_when_no_new_papers_are_requested(config, store, monkeypatch):
    rid, pipeline = setup_pipeline(config, store, monkeypatch, complete_rows(), [('known', paper())])
    store.finish_job(store.jobs(rid, 'paper')[0], [extracted()])
    def must_not_extract(*args):
        raise AssertionError('Cached paper was extracted again')
    monkeypatch.setattr(pipeline, 'paper', must_not_extract)
    _, counts = pipeline.execute(rid, '0')
    assert counts == {'measurements': 2, 'review': 0, 'all_observations': 3}
    pipeline.net.close()


def test_interruption_preserves_new_first_order_and_reuses_reference_lookups(config, store, monkeypatch):
    new = dict(paper('456'), doi='10.1234/new', title='Short title')
    rid, pipeline = setup_pipeline(config, store, monkeypatch, [record()], [('known', paper()), ('new', new)])
    calls = []
    def collect(payload, virus, run_id):
        calls.append(payload['pmcid'])
        if len(calls) == 2:
            raise KeyboardInterrupt()
        return [], []
    monkeypatch.setattr(pipeline, 'paper', collect)
    with pytest.raises(KeyboardInterrupt):
        pipeline.execute(rid, 'all')
    assert store.run(rid)['options']['selected_paper_keys'] == ['new', 'known']
    def no_new_lookup(*args):
        raise AssertionError('Reference metadata should be checkpointed')
    monkeypatch.setattr(pipeline.chembl, 'document', no_new_lookup)
    pipeline.execute(rid, None)
    assert calls == ['PMC456', 'PMC123', 'PMC123']
    pipeline.net.close()
