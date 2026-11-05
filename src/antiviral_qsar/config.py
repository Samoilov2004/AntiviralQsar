from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tomllib

PROMPT_VERSION = "2026-11-05-v1"
RULES_VERSION = "2026-11-05-v1"


def root_path() -> Path:
    return Path(os.environ.get("AQ_ROOT", Path(__file__).resolve().parents[2])).resolve()


def digest(value: object) -> str:
    text = value if isinstance(value, str) else json.dumps(value, sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(text.encode()).hexdigest()


class Config:
    def __init__(self, root: Path | None = None):
        self.root = (root or root_path()).resolve()
        self.settings = tomllib.loads((self.root / "config.toml").read_text())
        self.code_snapshot = {p.name:p.read_text() for p in sorted(Path(__file__).parent.glob('*.py'))}
        self.code_hash = digest(self.code_snapshot)
        for directory in ("data", "data/runs", "data/raw", "models", "logs", ".cache"):
            (self.root / directory).mkdir(parents=True, exist_ok=True)

    @property
    def model_url(self):
        model = self.settings["model"]
        return f"http://{model['host']}:{model['port']}"

    @property
    def model_dir(self):
        return self.root / "models" / "qwen3-8b-4bit"

    @property
    def model_revision(self):
        manifest = self.root / "models" / "manifest.json"
        return json.loads(manifest.read_text()).get("revision", "unknown") if manifest.exists() else "missing"

    @property
    def extractor_signature(self):
        return digest({"repo": self.settings["model"]["repo"], "revision": self.model_revision,
                       "prompt": PROMPT_VERSION, "rules": RULES_VERSION, "code_hash":self.code_hash,
                       "model_settings":self.settings['model']})

    def save_code_snapshot(self, rid):
        folder = self.root / 'data/runs' / rid / 'code' / self.code_hash[:16]
        folder.mkdir(parents=True, exist_ok=True)
        for name, content in self.code_snapshot.items():
            path = folder / name
            if not path.exists():
                path.write_text(content)
        return str(folder.relative_to(self.root))
