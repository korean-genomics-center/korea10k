#!/usr/bin/env bash
# Reproduce Supplementary Figure Sx end to end.
# Expects per_cpg_beta_by_genotype.with_values.filtered.tsv in the repo root,
# or run scripts/make_synthetic_input.py first to exercise the code path.
set -euo pipefail
PY="${PYTHON:-python3}"

echo "[1/5] parsing input"          && "$PY" scripts/01_parse.py
echo "[2/5] selecting sites"        && "$PY" scripts/00_select_sites.py
echo "[3/5] validating fast path"   && "$PY" scripts/02b_validate.py
echo "[4/5] simulating"             && "$PY" scripts/02_simulate.py --progress
echo "[5/5] building figure"        && "$PY" scripts/06_figure_power.py --extended

echo
echo "done:"
echo "  figures/fig_power_extended.pdf"
echo "  figures/fig_power_extended.png"
echo "  results/fig_power_source_data.tsv"
