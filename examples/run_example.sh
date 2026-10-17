#!/usr/bin/env bash
# Quick demo script to test bold-barcode-profiler on sample Culicidae COI contigs

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

echo "Running BOLD Barcode Profiler Example..."

python3 "${SCRIPT_DIR}/bold_barcode_profiler.py" \
    --query "${SCRIPT_DIR}/examples/example_contigs.fasta" \
    --bold "${SCRIPT_DIR}/examples/example_bold.tsv" \
    --outdir "${SCRIPT_DIR}/examples/example_results" \
    --prefix "demo_culicidae" \
    --marker "COI-5P,COII" \
    --min-identity 80 \
    --min-coverage 50 \
    --threads 2

echo "Example run finished! Results are in: ${SCRIPT_DIR}/examples/example_results"
