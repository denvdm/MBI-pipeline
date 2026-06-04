#!/usr/bin/env bash
set -euo pipefail

# Thin wrapper around the self-contained toy smoke test.
# Runs from the repo root regardless of where it is called from.
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

echo "[toy] repo root: $REPO_ROOT"
echo "[toy] running examples/toy/run_toy_example.R"
Rscript examples/toy/run_toy_example.R
