#!/usr/bin/env bash
# Recompute the compact workflow from included inputs.
#
# Default mode:
#   - recompute the structural and multispecies analyses;
#   - regenerate all diagnostic plots;
#   - validate archived expression derivatives from public datasets.
#
# Optional full expression mode downloads the accessioned large public files
# before rerunning the expression analysis:
#   bash scripts/run_all.sh --full-expression

set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python_bin="${RETROEVENT_PYTHON:-python3}"
mode="${1:-cached}"

if [[ "$mode" != "cached" && "$mode" != "--full-expression" ]]; then
    printf 'Usage: bash scripts/run_all.sh [--full-expression]\n' >&2
    exit 2
fi

"$python_bin" "$project_root/scripts/01_reconstruct_event.py" --project-root "$project_root"
"$python_bin" "$project_root/scripts/02_analyse_synteny.py"

if [[ "$mode" == "--full-expression" ]]; then
    bash "$project_root/scripts/fetch_expression_inputs.sh"
    "$python_bin" "$project_root/scripts/03_analyse_expression_full.py"
else
    "$python_bin" "$project_root/scripts/04_plot_expression_evidence.py"
fi

"$python_bin" "$project_root/scripts/03_validate_expression.py"
"$python_bin" "$project_root/scripts/05_validate_package.py"

printf 'All analyses and validation checks completed successfully.\n'
