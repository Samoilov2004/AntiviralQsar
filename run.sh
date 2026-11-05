#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
export AQ_ROOT="$ROOT"
export HF_HOME="$ROOT/.cache/huggingface"
export XDG_CACHE_HOME="$ROOT/.cache"
export MPLCONFIGDIR="$ROOT/.cache/matplotlib"
export PYTHONPATH="$ROOT/src${PYTHONPATH:+:$PYTHONPATH}"
export TOKENIZERS_PARALLELISM=false
export PYTHONUNBUFFERED=1
if [[ ! -x "$ROOT/.venv/bin/python" ]]; then
  echo "Run ./setup.sh first." >&2
  exit 1
fi
exec "$ROOT/.venv/bin/python" -m antiviral_qsar "$@"
