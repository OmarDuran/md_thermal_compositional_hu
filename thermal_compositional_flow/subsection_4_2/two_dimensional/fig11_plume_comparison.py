#!/usr/bin/env python
"""Q = 9 W/m^2 plume at one time: HU-PorePy TPFA, HU-PorePy MPFA and openCSMP.

Rows 1-2 are vapor saturation and fluid enthalpy maps, one column per scheme. Row 3 overlays the
Fig. 10A contours of all three (openCSMP solid with circles, TPFA dotted, MPFA dashed). Reads only; writes <out>.{png,pdf} here.

    python fig11_plume_comparison.py                    # 23.2 kyr
    python fig11_plume_comparison.py --years 25000 --csmp-run ../../opencsmp_data/two_dimensional/<tag>
"""
from __future__ import annotations

import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402

import fig10_gravity_comparison as F  # noqa: E402
import plot_csmp_snapshot as P  # noqa: E402

CSMP_RUN = "../../opencsmp_data/two_dimensional/hu_simplex_trunc_q9_dtmin1e-08"   # run_opencsmp.sh plume
MAP_XLIM = (-1.25, 1.25)       # the plume lies within about 1 km of the axis
FIELDS = ("pressure", "T_C", "s_l", "enthalpy")
P_LEVELS = [5.0, 10.0, 15.0, 25.0]    # fluid pressure [MPa]
H_LEVELS = [0.2, 1.0]                 # fluid enthalpy [MJ/kg]
H_LABELS = [(-1.62, 0.76), (-0.45, 1.25)]   # enthalpy label spots, clear of the pressure labels
TPFA_LS, MPFA_LS = (0, (1.2, 1.8)), (0, (5, 2.5))
MARK = dict(marker="o", ms=3.2, mfc="none", mew=0.8)   # openCSMP contour markers
MARK_SPACING = 0.15            # km of arc length between markers


def porepy_at(folder: str, years: float):
    """(Triangulation, nodal fields, cell fields) of the PorePy snapshot closest to ``years``."""
    t = np.array(json.load(open(os.path.join(F.HERE, folder, "times.json")))["time"]) / (365.0 * 86400.0)
    k = int(np.argmin(abs(t - years)))
    if abs(t[k] - years) > 1.0:
        raise SystemExit(f"{folder}: no snapshot at {years} yr (closest {t[k]:.0f} yr)")
    tri, cd = F.load(folder, k)
    return tri, {f: F.cell_to_point(tri, cd[f]) for f in FIELDS}, cd


def csmp_at(run: str, years: float):
    folder = os.path.join(F.HERE, run)
    snaps = dict(np.loadtxt(os.path.join(folder, "csmp_plume_snapshot_years.txt")).astype(int).tolist())
    idx = next((i for i, y in snaps.items() if abs(y - years) < 1.0), None)
    if idx is None:
        raise SystemExit(f"{run}: no snapshot at {years} yr")
    tri, f, _ = P.load_csmp(run, idx)
    return tri, f


def contours(ax, tri, f, ls, label=False, markers=False):
    """Fig. 10A contours of one scheme over the full width; small open circles in the line color."""
    kw = dict(linestyles=[ls], linewidths=1.3)
    sets = [ax.tricontour(tri, f["pressure"], levels=P_LEVELS, colors=F.P_COLOR, **kw),
            ax.tricontour(tri, f["T_C"], levels=F.T_LEVELS, colors=F.T_COLOR, **kw),
            ax.tricontour(tri, f["s_l"], levels=F.SL_LEVELS, colors=F.S_COLOR, **kw),
            ax.tricontour(tri, f["enthalpy"], levels=H_LEVELS, colors=F.H_COLOR, **kw)]
    if label:
        for cs, fmt, where in zip(sets, ("%g MPa", "%g°C", "%.2f", "%g MJ/kg"), (False, False, False, H_LABELS)):
            for tx in ax.clabel(cs, fmt=fmt, fontsize=9, inline=True, inline_spacing=5, manual=where):
                tx.set_bbox(dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    if markers:                     # after the labels, so the label gaps stay empty
        for cs, color in zip(sets, (F.P_COLOR, F.T_COLOR, F.S_COLOR, F.H_COLOR)):
            for path in cs.get_paths():
                for seg in _segments(path):
                    _circles(ax, seg, color)


def _segments(path):
    """Raw polylines of a path, split at its MOVETO codes (to_polygons simplifies them)."""
    v, c = path.vertices, path.codes
    if c is None:
        return [v]
    starts = np.flatnonzero(c == path.MOVETO)
    return np.split(v, starts[1:])


def _circles(ax, xy, color, spacing=MARK_SPACING):
    """Open circles every ``spacing`` km of arc length along one contour segment."""
    if len(xy) < 2:
        return
    s = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(xy, axis=0).T))])
    t = np.arange(0.5 * spacing, s[-1], spacing)
    ax.plot(np.interp(t, s, xy[:, 0]), np.interp(t, s, xy[:, 1]), ls="none", mec=color, **MARK, zorder=7)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--years", type=float, default=23200.0)
    ap.add_argument("--csmp-run", default=CSMP_RUN, help="openCSMP run folder (relative to here)")
    ap.add_argument("--out", default="fig11_plume_comparison")
    args = ap.parse_args()

    tri_t, pt_t, cd_t = porepy_at(F.TPFA_DIR, args.years)
    tri_m, pt_m, cd_m = porepy_at(F.MPFA_DIR, args.years)
    tri_c, f_c = csmp_at(args.csmp_run, args.years)

    cols = (("TPFA", "TPFA"), ("MPFA", "MPFA"), ("CSMP", "openCSMP"))
    sv = {"TPFA": cd_t["s_v"], "MPFA": cd_m["s_v"], "CSMP": f_c["s_v"]}
    hh = {"TPFA": cd_t["enthalpy"], "MPFA": cd_m["enthalpy"], "CSMP": f_c["enthalpy"]}
    tris = {"TPFA": tri_t, "MPFA": tri_m, "CSMP": tri_c}
    h_max = max(float(np.max(v)) for v in hh.values())

    fig = plt.figure(figsize=(9.0, 9.6), constrained_layout=True)
    axd = fig.subplot_mosaic([["sv_TPFA", "sv_MPFA", "sv_CSMP"],
                              ["h_TPFA", "h_MPFA", "h_CSMP"],
                              ["cont", "cont", "cont"]],
                             gridspec_kw={"height_ratios": [1.0, 1.0, 1.75]})
    maps = {}
    for j, (key, title) in enumerate(cols):
        for row, data, vmax in (("sv", sv, 1.0), ("h", hh, h_max)):
            ax = axd[f"{row}_{key}"]
            vals = data[key]
            if key == "CSMP":                      # nodal (vertex-centered) values
                pc = ax.tripcolor(tris[key], vals, cmap=F.SAT_CMAP, vmin=0.0, vmax=vmax,
                                  shading="gouraud", rasterized=True)
            else:                                  # cell values
                pc = ax.tripcolor(tris[key], facecolors=vals, cmap=F.SAT_CMAP, vmin=0.0, vmax=vmax,
                                  shading="flat", rasterized=True)
            ax.axvline(0.0, color="0.4", ls="--", lw=1.0, zorder=5)
            F._style_axes(ax, title if row == "sv" else None, ylabel=(j == 0), xlabel=(row == "h"))
            ax.set_xlim(*MAP_XLIM)
            ax.set_aspect("auto")
            if row == "h":
                ax.set_xlabel("distance [km]")
            maps[row] = pc
    fig.colorbar(maps["sv"], ax=[axd[f"sv_{k}"] for k, _ in cols], location="right",
                 fraction=0.03, pad=0.02, aspect=12).set_label("$s_v$")
    fig.colorbar(maps["h"], ax=[axd[f"h_{k}"] for k, _ in cols], location="right",
                 fraction=0.03, pad=0.02, aspect=12).set_label("$h$ [MJ/kg]")

    axc = axd["cont"]
    contours(axc, tri_t, pt_t, TPFA_LS)
    contours(axc, tri_m, pt_m, MPFA_LS)
    contours(axc, tri_c, f_c, "solid", label=True, markers=True)
    F._style_axes(axc, None)
    axc.set_aspect("auto")
    axc.set_title(f"{args.years / 1000:g} kyrs", loc="center")

    blank = Line2D([], [], color="none", label="")
    handles = [Line2D([], [], color=F.P_COLOR, lw=1.6, label="Pressure"),
               Line2D([], [], color="0.2", lw=1.6, label="openCSMP", mec="0.2", **MARK),
               Line2D([], [], color=F.T_COLOR, lw=1.6, label="Temperature"),
               Line2D([], [], color="0.2", lw=1.6, ls=TPFA_LS, label="HU-PorePy TPFA"),
               Line2D([], [], color=F.S_COLOR, lw=1.6, label="Liquid saturation"),
               Line2D([], [], color="0.2", lw=1.6, ls=MPFA_LS, label="HU-PorePy MPFA"),
               Line2D([], [], color=F.H_COLOR, lw=1.6, label="Enthalpy"), blank]
    axc.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=4,
               fontsize=12, frameon=True, columnspacing=1.4, handlelength=1.6)

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(F.HERE, f"{args.out}.{ext}"), dpi=300, bbox_inches="tight")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
