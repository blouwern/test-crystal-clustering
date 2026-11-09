#!/usr/bin/env bash
# Stage 2: generate the two multi-particle scenarios (signal + background)
# and evaluate the current reconstruction over the seed/wave/Nth grid.
set -euo pipefail
cd "$(dirname "$0")"
python3 generate.py
python3 evaluate.py
