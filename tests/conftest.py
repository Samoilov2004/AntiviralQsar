from pathlib import Path
import shutil
import pytest
from antiviral_qsar.config import Config
from antiviral_qsar.storage import Store


@pytest.fixture
def config(tmp_path):
    source = Path(__file__).resolve().parents[1] / "config.toml"
    shutil.copy(source, tmp_path / "config.toml")
    return Config(tmp_path)


@pytest.fixture
def store(config):
    return Store(config.root)
