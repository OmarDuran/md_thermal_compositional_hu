#!/usr/bin/env bash
# Section 4.1 of the paper: immiscible multiphase gravitational segregation.
#
# Usage:  ./run_subsection_4_1.sh [--list] [STEP ...]
#
# Steps (default: all, in this order):
#   checks    reduction-consistency table and monotonicity figure
#   column    column segregation with a passive phase (standalone solver, HU vs PPU)
#   barriers  barrier benchmark, fixed- and mixed-dimensional (standalone + PorePy)
#   collect   copy the paper figures to paper_figures/
#
# --list prints the commands without running them. Runs are strictly sequential (one
# process at a time). Finished runs are cached on disk and skipped by the barriers step;
# delete their output folders to recompute. Set PYTHON to choose the interpreter.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUB="${HERE}/subsection_4_1"
OUT="${HERE}/paper_figures"
PYTHON="${PYTHON:-python}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}" MPLBACKEND="${MPLBACKEND:-Agg}"

LIST=0
STEPS=()
for arg in "$@"; do
    case "$arg" in
        --list) LIST=1 ;;
        -h|--help) sed -n '2,15p' "$0"; exit 0 ;;
        checks|column|barriers|collect) STEPS+=("$arg") ;;
        *) echo "unknown step: $arg" >&2; exit 2 ;;
    esac
done
[ ${#STEPS[@]} -eq 0 ] && STEPS=(checks column barriers collect)

# Run a command inside subsection_4_1/, or only print it with --list.
run() {
    echo "+ (cd subsection_4_1 && $*)"
    [ "$LIST" -eq 1 ] && return 0
    (cd "$SUB" && "$@")
}

# Copy one produced file into paper_figures/, named as in the paper's images/ folder.
collect() {
    local src="$SUB/$1" dst="$OUT/$2"
    echo "+ cp subsection_4_1/$1 paper_figures/$2"
    [ "$LIST" -eq 1 ] && return 0
    if [ -f "$src" ]; then
        mkdir -p "$(dirname "$dst")" && cp "$src" "$dst"
    else
        echo "  missing: subsection_4_1/$1" >&2
    fi
}

step_checks() {
    # The table uses the coarse configuration (--quick); the monotonicity figure the full one.
    run "$PYTHON" completion_checks.py --quick
    run "$PYTHON" -c "import completion_checks as CC; CC.SAVE_PDF = True; CC.check_monotonicity_porepy()"
}

step_column() {
    # K = 10 D, ny = 200, 60 days (the script defaults).
    run "$PYTHON" column_segregation.py
}

step_barriers() {
    # Standalone PPU and HU-BM(mp) for N = 3, 4; PorePy HU fixed- and mixed-dimensional for
    # N = 3, 4; the equi-dimensional standalone runs for N = 3 last (the slow ones).
    run "$PYTHON" run_figures.py
}

step_collect() {
    local d=immiscible_multiphase_two_dimensional_flow
    collect figures/monotonicity.pdf "$d/monotonicity.pdf"
    collect tables/degenerate_equivalence.tex "$d/degenerate_equivalence.tex"
    collect passive_phase_example/column_segregation.pdf "$d/column_segregation_passive_phase.pdf"
    collect figures/n4/comparison_saturation_maps_fixed_dim_4_phases.pdf \
        "$d/comparison_saturation_maps_fixed_dim_4_phases.pdf"
    collect figures/n3/comparison_saturation_maps_mixed_dim_3_phases.pdf \
        "$d/comparison_saturation_maps_mixed_dim_3_phases.pdf"
}

for s in "${STEPS[@]}"; do
    echo "=== section 4.1: $s ($(date '+%Y-%m-%d %H:%M:%S'))"
    "step_$s"
done
