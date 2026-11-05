#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export AQ_ROOT="$ROOT"
export HF_HOME="$ROOT/.cache/huggingface"
export XDG_CACHE_HOME="$ROOT/.cache"
export PIP_CACHE_DIR="$ROOT/.cache/pip"
PYTHON="${AQ_PYTHON:-python3.11}"
if ! command -v "$PYTHON" >/dev/null; then
  echo "Python 3.11 is required. Set AQ_PYTHON=/path/to/python3.11 if needed." >&2
  exit 1
fi
mkdir -p "$ROOT/.cache" "$ROOT/models" "$ROOT/data" "$ROOT/logs"
[[ -x "$ROOT/.venv/bin/python" ]] || "$PYTHON" -m venv "$ROOT/.venv"
if [[ -f "$ROOT/requirements.lock.txt" ]]; then
  "$ROOT/.venv/bin/python" -m pip install -r "$ROOT/requirements.lock.txt"
else
  "$ROOT/.venv/bin/python" -m pip install "$ROOT[mlx,dev]"
fi
"$ROOT/.venv/bin/python" -m pip install --no-deps -e "$ROOT"
"$ROOT/.venv/bin/python" -m pip freeze --exclude-editable > "$ROOT/requirements.lock.txt"
"$ROOT/run.sh" model download
"$ROOT/run.sh" doctor
