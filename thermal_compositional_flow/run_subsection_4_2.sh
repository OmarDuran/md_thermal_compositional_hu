#!/usr/bin/env bash
# Sections 4.2 and 4.3 of the paper and the figures of Appendix B, whose scripts all live in
# subsection_4_2/.
#
# Usage:  ./run_subsection_4_2.sh [--list] [STEP ...]
#
# Steps (default: all, in this order):
#   weis1d    one-dimensional Weis et al. (2014) benchmarks: Fig. 5 (Section 4.2) and the
#             appendix Figs. 4 and 6, standalone solver plus the PorePy overlays
#   plume     two-dimensional geothermal plume, TPFA vs MPFA gravity (Section 4.2)
#   co2       water-CO2 flow across the fracture network, three runs (Section 4.3)
#   phz       water-CO2 PHZ diagrams of the OBL tables (Appendix B)
#   collect   copy the paper figures to paper_figures/
#
# --list prints the commands without running them. Runs are strictly sequential (one
# process at a time). The weis1d step reuses its _cache/ (delete it to recompute); the
# plume and co2 steps always rerun. Set PYTHON to choose the interpreter.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUB="${HERE}/subsection_4_2"
OUT="${HERE}/paper_figures"
PYTHON="${PYTHON:-python}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}" MPLBACKEND="${MPLBACKEND:-Agg}"

LIST=0
STEPS=()
for arg in "$@"; do
    case "$arg" in
        --list) LIST=1 ;;
        -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
        weis1d|plume|co2|phz|collect) STEPS+=("$arg") ;;
        *) echo "unknown step: $arg" >&2; exit 2 ;;
    esac
done
[ ${#STEPS[@]} -eq 0 ] && STEPS=(weis1d plume co2 phz collect)

# Run a command inside subsection_4_2/<dir>, or only print it with --list.
run() {
    local dir="$1"; shift
    echo "+ (cd subsection_4_2/$dir && $*)"
    [ "$LIST" -eq 1 ] && return 0
    (cd "$SUB/$dir" && "$@")
}

# Copy one produced file into paper_figures/, named as in the paper's images/ folder.
collect() {
    local src="$SUB/$1" dst="$OUT/$2"
    echo "+ cp subsection_4_2/$1 paper_figures/$2"
    [ "$LIST" -eq 1 ] && return 0
    if [ -f "$src" ]; then
        mkdir -p "$(dirname "$dst")" && cp "$src" "$dst"
    else
        echo "  missing: subsection_4_2/$1" >&2
    fi
}

step_weis1d() {
    # --serial keeps the 1-D sweeps in one process (a parallel pool loads the tables per worker).
    run one_dimensional "$PYTHON" run_workflow.py --serial --pdf
}

step_plume() {
    # Truncated domain, simplex mesh (100 m), Q = 9 W/m^2, 25 kyr: TPFA, then MPFA.
    run two_dimensional "$PYTHON" porepy_2d_solver.py --grid-type simplex --truncated-domain --q-anomaly 9
    run two_dimensional "$PYTHON" porepy_2d_solver.py --grid-type simplex --truncated-domain --q-anomaly 9 \
        --consistent
    run two_dimensional "$PYTHON" fig10_gravity_comparison.py
}

step_co2() {
    # TPFA, MPFA, and MPFA without gravity, 50 days each, then the two figures.
    run two_dimensional "$PYTHON" run_co2_md.py
}

step_phz() {
    run two_dimensional/h2o_co2_obl "$PYTHON" co2_phz_diagrams.py
}

step_collect() {
    collect one_dimensional/figures/reproduce_weis_fig_5.pdf weis_one_dimensional_flow/reproduce_weis_fig_5.pdf
    collect one_dimensional/figures/reproduce_weis_fig_4.pdf weis_one_dimensional_flow/reproduce_weis_fig_4.pdf
    collect one_dimensional/figures/reproduce_weis_fig_6.pdf weis_one_dimensional_flow/reproduce_weis_fig_6.pdf
    collect two_dimensional/fig10_gravity_comparison.pdf weis_two_dimensional_flow/fig10_gravity_comparison.pdf
    collect two_dimensional/co2_rho_h_contours.pdf co2_two_dimensional_flow/co2_rho_h_contours.pdf
    collect two_dimensional/co2_breakthrough.pdf co2_two_dimensional_flow/co2_breakthrough.pdf
    collect two_dimensional/h2o_co2_obl/co2_phz_diagrams.pdf co2_two_dimensional_flow/co2_phz_diagrams.pdf
}

for s in "${STEPS[@]}"; do
    echo "=== subsection_4_2: $s ($(date '+%Y-%m-%d %H:%M:%S'))"
    "step_$s"
done
