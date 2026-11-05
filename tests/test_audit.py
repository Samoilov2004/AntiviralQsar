import csv
import json
from typer.testing import CliRunner
from antiviral_qsar import cli
from antiviral_qsar.chemistry import validate_measurement
from antiviral_qsar.pipeline import export_run
from test_integrity import measurement


def test_blank_audit_cannot_validate_extractor(config, store, monkeypatch):
    monkeypatch.setattr(cli, 'components', lambda:(config,store))
    rid=store.new_run({'name':'Dengue virus'}, {})
    store.add_job(rid, 'europepmc','paper','p',{})
    m=validate_measurement(measurement(extraction_method='llm',extractor_signature='fixture'))
    store.finish_job(store.jobs(rid)[0], [m])
    runner=CliRunner()
    assert runner.invoke(cli.app,['audit','prepare',rid,'--size','100']).exit_code == 0
    path=config.root/'data/runs'/rid/'audit-llm-42.csv'
    result=runner.invoke(cli.app,['audit','evaluate',str(path)])
    assert result.exit_code != 0
    assert not store.quality_passed('fixture')
    with path.open(encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==1 and rows[0]['correct']==''
    assert json.loads(path.with_suffix('.json').read_text())['signature']=='fixture'


def test_revoked_quality_gate_restores_review_status(config, store):
    rid=store.new_run({'name':'Dengue virus'}, {})
    store.add_job(rid,'europepmc','paper','p',{})
    # A previously approved record remains stored without a review reason.
    m=validate_measurement(measurement(extraction_method='llm',extractor_signature='fixture'),llm_gate=True)
    store.finish_job(store.jobs(rid)[0],[m])
    store.db.execute('INSERT INTO quality_gates VALUES(?,?,?)',('fixture',json.dumps({'passed':True}),'test'))
    store.db.commit()
    assert export_run(config,store,rid)[1]['measurements']==1
    store.db.execute('UPDATE quality_gates SET data=?', (json.dumps({'passed':False}),))
    store.db.commit()
    counts=export_run(config,store,rid)[1]
    assert counts['measurements']==0 and counts['review']==1
