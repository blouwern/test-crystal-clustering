#!/usr/bin/env bash
# Stage 3: generate separation-binned 2gamma pairs and measure merge rate vs distance.
set -euo pipefail
cd "$(dirname "$0")"
python3 generate.py
python3 evaluate.py
