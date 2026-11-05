import json
from unittest.mock import Mock
import pytest
import huggingface_hub
from antiviral_qsar.model import LocalModel


SHA = '545dc4251c05440727734bcd94334791f6ab0192'


def test_fresh_install_downloads_exact_committed_model(config, store, monkeypatch):
    repo = config.settings['model']['repo']
    (config.root/'model.lock.json').write_text(json.dumps({'repo':repo,'revision':SHA}))
    download = Mock()
    monkeypatch.setattr(huggingface_hub, 'snapshot_download', download)
    assert LocalModel(config, store).download() == SHA
    assert download.call_args.kwargs['revision'] == SHA
    assert download.call_args.kwargs['local_dir'] == config.model_dir
    assert config.model_revision == SHA


def test_different_installed_model_is_not_overwritten(config, store, monkeypatch):
    repo = config.settings['model']['repo']
    (config.root/'model.lock.json').write_text(json.dumps({'repo':repo,'revision':SHA}))
    installed = {'repo':repo,'revision':'a'*40}
    (config.root/'models/manifest.json').write_text(json.dumps(installed))
    download = Mock()
    monkeypatch.setattr(huggingface_hub, 'snapshot_download', download)
    with pytest.raises(ValueError, match='Installed model differs'):
        LocalModel(config, store).download()
    download.assert_not_called()
    assert config.model_revision == 'a'*40
