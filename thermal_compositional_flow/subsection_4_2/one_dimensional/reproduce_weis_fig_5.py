"""Weis (2014) Figure 5 -- two-phase, pure water. 2x2 panels: columns {horizontal 200 yr,
vertical 1000 yr}, rows {temperature+pressure, liquid saturation}. The single weis brine engine at
z=0 (``weis_1d_solver.run_brine(**FIG5)``) for PPU / HU / HU-mwp, overlaid on the openCSMP
reference (csmp_1d_solver.py). The heavy N=800 runs are cached in ``_cache/`` (delete to recompute).

    python reproduce_weis_fig_5.py            # compute (or load cache) + render figures/reproduce_weis_fig_5
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_weis_common as C  # noqa: E402  (pins BLAS threads on import)
import weis_1d_solver as m   # noqa: E402
import plot_style as ps      # noqa: E402

N = 800
TAG = "fig5"
CASES = ("horizontal", "vertical")
YEARS = {"horizontal": "200", "vertical": "1000"}

# openCSMP on the Weis (2014) grid of N=200 elements, next to the openCSMP reference (N=800).
CSMP_COARSE_LABEL = r"openCSMP ($\Delta x = 10$ m)"

# PorePy approximation overlay: the converged HU profile from porepy_1d_solver, cached as
# _cache/fig5_porepy_{case}_hu_N800_l<level>.pkl. Drawn as black x markers over the weis-HU reference so the
# agreement (or drift) reads directly.
POREPY_C = "black"
POREPY_LABEL = r"HU-PorePy"
# If a porepy overlay pickle is missing, run porepy_1d_solver to generate it (so one `python
# reproduce_weis_fig_5.py` produces the porepy data AND the figure). HEAVY -- the PorePy solve is minutes
# per case (vertical is 1000 yr). Set False to only plot porepy caches that already exist.
AUTORUN_POREPY = True


def _load_porepy(case, scheme="hu", N=N):
    """porepy_1d_solver profile pickle (y[m], T[K], p[MPa], s_liq) for ``case``/``scheme``, normalised
    to the SI plot convention (p -> Pa) that ``ps.to_plot_units`` consumes. If the pickle is missing
    and ``AUTORUN_POREPY``, run ``porepy_1d_solver.run_case`` to make it (lazy import, so a warm-cache
    re-plot never imports porepy). Returns the dict, or None if unavailable."""
    import pickle
    path = os.path.join(C.CACHE_DIR, f"fig5_porepy_{case}_{scheme}_N{N}_l{m.TABLE_LEVEL}.pkl")
    if not os.path.exists(path) and AUTORUN_POREPY:
        try:
            import porepy_1d_solver as pp1d                 # lazy: only imports porepy on a cold cache
            if N == pp1d.N_CELLS:
                print(f"[fig5] porepy overlay cache missing for {case}/{scheme} -- running "
                      f"porepy_1d_solver.run_case (heavy: PorePy solve) ...", flush=True)
                C.run_in_own_process("import porepy_1d_solver as P; P.run_case("
                                     f"{case!r}, weighted_perm={scheme == 'hu_mwp'})")  # writes the same pickle path
            else:
                print(f"[fig5] porepy overlay skipped for {case}/{scheme}: fig N={N} "
                      f"!= porepy_1d_solver N={pp1d.N_CELLS}", flush=True)
        except Exception as exc:                            # never let an overlay break the figure
            print(f"[fig5] porepy overlay auto-run failed for {case}/{scheme}: {exc}", flush=True)
    if not os.path.exists(path):
        return None
    with open(path, "rb") as f:
        d = dict(pickle.load(f))
    d["p"] = d["p"] * 1.0e6                                # MPa (PorePy native) -> Pa (SI, as weis)
    print(f"[fig5-porepy] ({scheme!r}, {case!r})   cached total_it={int(d['total_it'])}", flush=True)
    return d


def compute(N=N, parallel=True, skip=frozenset()):
    out = C.sweep(TAG, list(CASES), m.FIG5, N, parallel=parallel,            # PPU / HU / HU-mwp (minus skipped)
                  schemes=C.active_schemes(skip))
    return out


def plot(out, stem="reproduce_weis_fig_5", skip=frozenset()):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    ps.apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(ps.TEXTWIDTH_IN, 5.6), sharex="col")
    tags = (("(a)", "(b)"), ("(c)", "(d)"))
    panels = []
    has_pp = False                                       # True once a PorePy overlay is actually drawn
    for j, case in enumerate(CASES):
        ax_tp, ax_s = axes[0, j], axes[1, j]
        ax_p = ax_tp.twinx(); ax_p.grid(False)
        res = {sk: out[(sk, case)] for sk in ps.SCHEMES if (sk, case) in out}
        panels.append(res)
        C.draw_tp(ax_tp, ax_p, res,
                  ref_T=C.csmp_ref(f"fig_5_{case}", "T"),
                  ref_p=C.csmp_ref(f"fig_5_{case}", "p"), ref_style=C.CSMP_STYLE)
        C.draw_s(ax_s, res, ref_s=C.csmp_ref(f"fig_5_{case}", "s_liq"), ref_style=C.CSMP_STYLE)
        extra_it = []
        for ax, fld, col in ((ax_tp, "T", ps.CSMP_COARSE["warm"]), (ax_p, "p", ps.CSMP_COARSE["cool"]),
                             (ax_s, "s_liq", ps.CSMP_COARSE["warm"])):   # openCSMP N=200 (no Newton count)
            ref = C.csmp_ref(f"fig_5_{case}_n200", fld)
            if ref is not None:
                ax.plot(*ref, color=col, ls=ps.CSMP_COARSE["dash"], lw=ps.CURVE_LW, zorder=3)
        pp_res = None if C.is_skipped("hu-porepy", skip) else _load_porepy(case)  # PorePy HU overlay
        if pp_res is not None:                                # warm T / cool p / warm s_liq, PorePy dash
            has_pp = True
            for ax, fld, col in ((ax_tp, "T", ps.POREPY["warm"]), (ax_p, "p", ps.POREPY["cool"]),
                                 (ax_s, "s_liq", ps.POREPY["warm"])):
                x, y = ps.to_plot_units(pp_res, fld)
                ax.plot(x, y, color=col, ls=ps.POREPY["dash"], lw=ps.CURVE_LW, zorder=6)
            extra_it.append((ps.POREPY["warm"], ps.POREPY["dash"], pp_res["total_it"]))
        C.iteration_legend(ax_s, res, loc="lower right", extra=extra_it)
        ax_tp.set_title(fr"\textbf{{{case.capitalize()} (${YEARS[case]}$ yr)}}")
        ax_tp.set_xlim(0.0, 2.0)
        ps.panel_tag(ax_tp, tags[0][j], loc=(0.04, 0.09), va="bottom")   # T+p high at top-left -> tag low
        ps.panel_tag(ax_s, tags[1][j])                                    # s_liq low at top-left -> tag high
        # left column carries the T / s_liq axes; right column carries the pressure axis
        if j == 0:
            ax_tp.set_ylabel(ps.FIELD_LABEL["T"], color=C.WEIS_T)
            ax_tp.tick_params(axis="y", colors=C.WEIS_T)
            ax_s.set_ylabel(ps.FIELD_LABEL["s_liq"], color=C.WEIS_T)   # liquid sat = warm (left)
            ax_s.tick_params(axis="y", colors=C.WEIS_T)
            ax_p.tick_params(axis="y", labelright=False)
        else:
            ax_p.set_ylabel(ps.FIELD_LABEL["p"], color=C.WEIS_P)
            ax_p.tick_params(axis="y", colors=C.WEIS_P)
            ax_tp.tick_params(axis="y", labelleft=False)
            ax_s.tick_params(axis="y", labelleft=False)
        ax_s.set_xlabel(ps.DIST_LABEL)

    # single bottom legend: scheme colours (names only) + the T/p and reference key. Per-case
    # iteration counts live in each panel's small legend (they differ by orientation).
    handles = C.scheme_handles(only=C.active_schemes(skip))
    handles.append(Line2D([0], [0], color=ps.CSMP_COARSE["warm"], ls=ps.CSMP_COARSE["dash"], lw=ps.CURVE_LW, label=CSMP_COARSE_LABEL))
    if has_pp:                                          # key the overlay only when it was actually drawn
        handles.append(Line2D([0], [0], color=ps.POREPY["warm"], ls=ps.POREPY["dash"], lw=ps.CURVE_LW, label=POREPY_LABEL))
    handles += [Line2D([0], [0], color=ps.SCHEMES["hu"]["warm"], ls="-", lw=3, label=r"$T$ (warm, left)"),
                Line2D([0], [0], color=ps.SCHEMES["hu"]["cool"], ls="-", lw=3, label=r"$p$ (cool, right)"),
                C.ref_legend_handle(C.CSMP_LABEL, C.CSMP_STYLE)]
    fig.tight_layout()
    ps.bottom_legend(fig, handles, [h.get_label() for h in handles], ncol=4)
    ps.savefig(fig, stem, C.OUT_DIR)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Weis (2014) Fig 5 (two-phase pure water).")
    ap.add_argument("--skip-solver", default="", dest="skip",
                    help="comma-separated solvers to skip: ppu, hu, hu-mwp, hu-porepy")
    args = ap.parse_args(argv)
    skip = C.parse_skip(args.skip)
    plot(compute(skip=skip), skip=skip)


if __name__ == "__main__":
    main()
