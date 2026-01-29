#!/usr/bin/env bash
set -euo pipefail

# Run from repo root regardless of where called from
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

CFG="examples/toy/config_toy.yaml"

echo "[toy] Using config: $CFG"
echo "[toy] Repo root: $REPO_ROOT"

# Option A: if you have a single orchestrator script
if [[ -f scripts/run_all.sh ]]; then
  echo "[toy] Running scripts/run_all.sh (may need to adjust flags)..."
  bash scripts/run_all.sh --config "$CFG"
  echo "[toy] Done."
  exit 0
fi

# Option B: fall back to showing available scripts
echo "[toy] No scripts/run_all.sh found. Available scripts:"
ls -la scripts || true
echo "[toy] Please edit examples/toy/run_toy.sh to call your preferred entry point."
exit 1
