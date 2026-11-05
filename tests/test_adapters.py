import pytest
from antiviral_qsar.sources import PubChem, pubchem_unit
from antiviral_qsar.pipeline import deduplicate
from test_integrity import measurement
from antiviral_qsar.chemistry import validate_measurement, normalize


class FakeNetwork:
    def __init__(self, payload):
        self.payload = payload

    def json(self, url, params=None):
        if url.endswith('/sids/JSON'):
            return {'InformationList': {'Information': [{'SID': [100, 101]}]}}, 'sids.json'
        return self.payload, 'results.json'


@pytest.mark.parametrize('wrapper', ['PC_AssaySubmit', 'PC_AssayContainer'])
def test_pubchem_single_and_multiple_submission_shapes(wrapper):
    submission = {'data': [{'sid': 100, 'data': []}, {'sid': 101, 'data': []}]}
    value = submission if wrapper == 'PC_AssaySubmit' else [submission]
    assert len(list(PubChem(FakeNetwork({wrapper: value})).results('1'))) == 2


def test_pubchem_malformed_results_never_silently_drop_measurements():
    with pytest.raises(ValueError, match='unrecognized'):
        list(PubChem(FakeNetwork({'unexpected': []})).results('1'))
    with pytest.raises(ValueError, match='incomplete'):
        list(PubChem(FakeNetwork({'PC_AssaySubmit': {'data': [{'sid': 100}]}})).results('1'))


def test_mass_unit_is_not_molar():
    assert pubchem_unit(9) == 'mg/mL'
    assert pubchem_unit(14) == 'M'
    with pytest.raises(ValueError, match='unsupported_units'):
        normalize('2', pubchem_unit(9))


def test_equal_replicates_within_source_are_retained():
    first = validate_measurement(measurement(source='chembl'))
    second = validate_measurement(measurement(source='chembl', observation_id='2', source_record_id='2'))
    assert not any(m.duplicate_of for m in deduplicate([first, second]))


def test_prefetch_keeps_sid_to_cid_mapping():
    class Net:
        def json(self, url, params=None):
            if '/substance/' in url:
                return {'InformationList': {'Information': [{'SID': 1, 'CID': [99]}, {'SID': 2, 'CID': [20, 30]}]}}, 'raw'
            if '/cid/99/' in url:
                return {'PropertyTable': {'Properties': [{'CID': 99, 'SMILES': 'CCO'}]}}, 'raw'
            raise AssertionError(url)
    source = PubChem(Net())
    source.prefetch_substances([1,2])
    assert source.sid_compound(1)[0]['CID'] == 99
    assert source.sid_compound(2)[1] == 'ambiguous_substance_structure'


def test_later_table_chunks_do_not_replay_data_rows_as_headers():
    from antiviral_qsar.documents import row_chunks
    rows = [['Compound', 'EC50']] + [[f'compound{i}', i] for i in range(12)]
    chunks = list(row_chunks(rows, 't'))
    assert len(chunks) == 2
    assert 'compound0' not in chunks[1].text


def test_pipeline_resume_after_interrupt_does_not_repeat_completed_assay(config, store, monkeypatch):
    from antiviral_qsar.pipeline import Pipeline
    from antiviral_qsar.schema import Virus
    virus = Virus(query='Dengue virus', name='Dengue virus', tax_id='12637')
    rid = store.new_run(virus.model_dump(), {'assay_limit': 2})
    store.add_jobs(rid, 'chembl', 'assay', [('A', {'assay_chembl_id':'A'}), ('B', {'assay_chembl_id':'B'})])
    pipeline = Pipeline(config, store, use_model=False)
    calls = []
    def assay(payload, virus):
        calls.append(payload['assay_chembl_id'])
        if calls == ['A', 'B']:
            raise KeyboardInterrupt()
        return [validate_measurement(measurement(observation_id=payload['assay_chembl_id']))]
    monkeypatch.setattr(pipeline, 'chembl_assay', assay)
    with pytest.raises(KeyboardInterrupt):
        pipeline.execute(rid, '0')
    assert [j['status'] for j in store.jobs(rid)] == ['done', 'running']
    pipeline.execute(rid, '0')
    assert calls == ['A', 'B', 'B']
    assert len(store.observations(rid)) == 2
    pipeline.net.close()


def test_host_organism_is_not_exported_as_virus(config, store):
    from antiviral_qsar.pipeline import Pipeline
    from antiviral_qsar.schema import Virus
    pipeline = Pipeline(config, store, use_model=False)
    class Chembl:
        def activities(self, aid):
            return [({'activity_id':1, 'molecule_chembl_id':'M1', 'standard_type':'EC50',
                      'type':'EC50', 'value':'1', 'units':'uM'}, 'raw')]
        def document(self, did):
            return {}
        def molecule(self, cid):
            return {'molecule_type':'Small molecule', 'molecule_structures':{'canonical_smiles':'CCO'}}
        def prefetch_molecules(self, ids):
            pass
    pipeline.chembl = Chembl()
    virus = Virus(query='Dengue virus',name='Dengue virus',tax_id='12637')
    data = pipeline.chembl_assay({'assay_chembl_id':'A', 'assay_organism':'Homo sapiens',
            'assay_tax_id':9606, 'assay_cell_type':'Huh7',
            'description':'Antiviral activity against Dengue virus in infected Huh7 cells'}, virus)
    assert data[0].virus == 'Dengue virus'
    assert data[0].virus_tax_id == '12637'
    pipeline.net.close()


def test_resume_preserves_original_paper_selection(config, store, monkeypatch):
    from antiviral_qsar.pipeline import Pipeline
    from antiviral_qsar.schema import Virus
    virus = Virus(query='Dengue virus', name='Dengue virus', tax_id='12637')
    rid = store.new_run(virus.model_dump(), {'papers':'2'})
    store.add_jobs(rid,'europepmc','paper',[(str(i), {'pmcid':f'PMC{i}', 'isOpenAccess':'Y'}) for i in range(3)])
    pipeline = Pipeline(config, store, use_model=False)
    calls = []
    def paper(payload, virus, rid):
        calls.append(payload['pmcid'])
        if len(calls) == 2:
            raise KeyboardInterrupt()
        return [], []
    monkeypatch.setattr(pipeline, 'paper', paper)
    with pytest.raises(KeyboardInterrupt):
        pipeline.execute(rid, '2')
    pipeline.execute(rid, None)
    assert calls == ['PMC0', 'PMC1', 'PMC1']
    assert [j['status'] for j in store.jobs(rid)] == ['done', 'done', 'pending']
    pipeline.net.close()


def test_pdf_fallback_uses_only_declared_open_pdf():
    from antiviral_qsar.sources import EuropePMC
    class Net:
        def get(self, url, **kwargs):
            if url.endswith('/fullTextXML'):
                raise RuntimeError('HTTP_404')
            assert url == 'https://example.org/open.pdf'
            return b'%PDF-1.7 example', 'pdf-snapshot'
    paper = {'pmcid':'PMC1','isOpenAccess':'Y','fullTextUrlList':{'fullTextUrl':[
        {'documentStyle':'pdf','availability':'Subscription','url':'https://example.org/paywall.pdf'},
        {'documentStyle':'pdf','availability':'Open access','url':'https://example.org/open.pdf'}]}}
    assert EuropePMC(Net()).fulltext(paper)[1] == 'pdf-snapshot'
