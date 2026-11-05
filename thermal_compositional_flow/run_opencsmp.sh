#!/usr/bin/env bash
# openCSMP (CSMP++) reference runs of the paper: the one-dimensional Weis et al. (2014)
# benchmarks of Section 4.2 and Appendix B, and the two-dimensional plume of Section 4.2.
# Run it before run_subsection_4_2.sh, whose figures read these profiles.
#
# Usage:  ./run_opencsmp.sh [--list] [STEP ...]
#
# All output goes to opencsmp_data/, which the plotting scripts read:
#   opencsmp_data/one_dimensional/<case>/csmp_profile_1.csv    fig_weis_common.csmp_ref
#   opencsmp_data/two_dimensional/<tag>/csmp_plume_*.vtu         2D plume snapshots
#
# Steps (default: all, in this order):
#   weis1d    1D columns on dx = 2.5 m (Fig. 4 cases, Fig. 5) and dx = 10 m (Fig. 5 coarse,
#             Fig. 6 water and water-salt)
#   plume     2D plume, simplex mesh (100 m), Q = 9 W/m^2, 25 kyr
#   save      pack opencsmp_data/ as opencsmp_data.tar.gz (without the lookup-table links
#             and the long scheme logs)
#
# Optional step, not run by default (no figure reads it):
#   extras    grid checks, Fig. 4 on dx = 10 m, Fig. 6 on dx = 2.5 m (the water-salt case
#             collapses near 215 yr, as reported in Appendix B), and the plume on a 72 m mesh
#             (about as many vertices as the 100 m mesh has cells)
#
# opencsmp_data/ is ignored by git; this script regenerates all of it. --list prints the
# commands without running them. Runs are strictly sequential. Set PYTHON to choose the
# interpreter.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUB="${HERE}/subsection_4_2"
OUT="${HERE}/opencsmp_data"
PYTHON="${PYTHON:-python}"

LIST=0
STEPS=()
for arg in "$@"; do
    case "$arg" in
        --list) LIST=1 ;;
        -h|--help) sed -n '2,28p' "$0"; exit 0 ;;
        weis1d|plume|extras|save) STEPS+=("$arg") ;;
        *) echo "unknown step: $arg" >&2; exit 2 ;;
    esac
done
[ ${#STEPS[@]} -eq 0 ] && STEPS=(weis1d plume save)

# Run a command inside subsection_4_2/<dir>, or only print it with --list.
run() {
    local dir="$1"; shift
    echo "+ (cd subsection_4_2/$dir && $*)"
    [ "$LIST" -eq 1 ] && return 0
    (cd "$SUB/$dir" && "$@")
}

step_weis1d() {
    # dx = 2.5 m (800 elements). Fig. 6 runs on dx = 10 m: on 2.5 m the water-salt case
    # collapses near 215 yr (csmp_1d_solver.py sets 200 elements for it).
    run one_dimensional "$PYTHON" csmp_1d_solver.py \
        fig_4_hp_horizontal fig_4_hp_vertical fig_4_mp_horizontal fig_4_mp_vertical \
        fig_4_lp_horizontal fig_4_lp_vertical fig_5_horizontal fig_5_vertical
    run one_dimensional "$PYTHON" csmp_1d_solver.py --elements 200 \
        fig_5_horizontal fig_5_vertical fig_6_pw fig_6_salt
}

step_plume() {
    # dt_min 1e-8 yr lets the scheme recover from transient CFL cuts.
    run two_dimensional "$PYTHON" csmp_2d_solver.py --grid-type simplex --truncated-domain \
        --q-anomaly 9 --dt-min 1e-8
}

step_extras() {
    run one_dimensional "$PYTHON" csmp_1d_solver.py --elements 200 \
        fig_4_hp_horizontal fig_4_hp_vertical fig_4_mp_horizontal fig_4_mp_vertical \
        fig_4_lp_horizontal fig_4_lp_vertical
    run one_dimensional "$PYTHON" csmp_1d_solver.py fig_6_pw
    run one_dimensional "$PYTHON" csmp_1d_solver.py --elements 800 --cfl-scaling 0.01 fig_6_salt \
        || echo "fig_6_salt on 800 elements stopped, as expected"
    run two_dimensional "$PYTHON" csmp_2d_solver.py --grid-type simplex --truncated-domain \
        --q-anomaly 9 --dt-min 1e-8 --cell-size 72
}

step_save() {
    echo "+ tar -czf opencsmp_data.tar.gz opencsmp_data"
    [ "$LIST" -eq 1 ] && return 0
    cat /workdir/openCSMP/OPENCSMP_VERSION > "$OUT/OPENCSMP_VERSION" 2>/dev/null || true
    tar -czf "${OUT}.tar.gz" -C "$HERE" --exclude='*.bin' --exclude='csmp.log' \
        --exclude='csmp_plume.log' --exclude='Runtime.log' "$(basename "$OUT")"
    du -sh "$OUT" "${OUT}.tar.gz"
}

for s in "${STEPS[@]}"; do
    echo "=== openCSMP: $s ($(date '+%Y-%m-%d %H:%M:%S'))"
    "step_$s"
done
