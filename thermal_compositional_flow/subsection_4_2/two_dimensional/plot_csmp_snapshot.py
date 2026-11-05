#!/usr/bin/env python
"""One openCSMP plume snapshot in the style of fig10_gravity_comparison.py.

Vapor saturation and fluid enthalpy maps, and the Fig. 10A contours with the final PorePy MPFA
snapshot dashed for reference. Reads only; writes <out>.{png,pdf} here.

    python plot_csmp_snapshot.py --index 464
"""
from __future__ import annotations

import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import meshio  # noqa: E402
import numpy as np  # noqa: E402

import fig10_gravity_comparison as F  # noqa: E402


def load_csmp(run: str, index: int):
    """(Triangulation [km], nodal fields in plot units, years) of one openCSMP snapshot."""
    folder = os.path.join(F.HERE, run)
    m = meshio.read(os.path.join(folder, f"csmp_plume_Model_(REGION)_{index}.vtu"))
    tri = mtri.Triangulation((m.points[:, 0] - F.X_PLUME) / 1000.0,
                             (F.Y_SURFACE - m.points[:, 1]) / 1000.0,
                             triangles=m.cells_dict["triangle"])
    pd = {k: np.asarray(v, float).ravel() for k, v in m.point_data.items()}
    f = {"pressure": pd["fluid pressure"] / 1e6, "T_C": pd["temperature"],
         "s_l": pd["saturation liquid"], "s_v": pd["saturation vapor"],
         "enthalpy": pd["fluid enthalpy"] / 1e6}
    years = dict(np.loadtxt(os.path.join(folder, "csmp_plume_snapshot_years.txt")).astype(int).tolist())
    return tri, f, years.get(index)


def contours(ax, tri, f, ls, label):
    kw = dict(linestyles=ls, linewidths=1.3)
    out = [ax.tricontour(tri, f["pressure"], levels=F.P_LEVELS, colors=F.P_COLOR, **kw),
           ax.tricontour(tri, f["T_C"], levels=F.T_LEVELS, colors=F.T_COLOR, **kw),
           ax.tricontour(tri, f["s_l"], levels=F.SL_LEVELS, colors=F.S_COLOR, **kw),
           ax.tricontour(tri, f["enthalpy"], levels=F.H_LEVELS, colors=F.H_COLOR, **kw)]
    if label:
        for cs, fmt in zip(out, ("%g MPa", "%g°C", "%.2f", "%.1f MJ/kg")):
            F._clabel(ax, cs, fmt)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default="../../opencsmp_data/two_dimensional/hu_simplex_trunc_q9_dtmin1e-08")
    ap.add_argument("--index", type=int, default=464)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    tri, f, yrs = load_csmp(args.run, args.index)
    tri_m, cd_m = F.load(F.MPFA_DIR, None)
    pp = {k: F.cell_to_point(tri_m, cd_m[k]) for k in ("pressure", "T_C", "s_l", "enthalpy")}

    fig = plt.figure(figsize=(9.0, 9.0), constrained_layout=True)
    axd = fig.subplot_mosaic([["sv", "h"], ["cont", "cont"]], gridspec_kw={"height_ratios": [1.0, 1.9]})
    for key, fld, lab, vmax in (("sv", "s_v", "vapor saturation $s_v$ [-]", 1.0),
                                ("h", "enthalpy", "fluid enthalpy $h$ [MJ/kg]", None)):
        ax = axd[key]
        tpc = ax.tripcolor(tri, f[fld], cmap=F.SAT_CMAP, vmin=0.0, vmax=vmax,
                           shading="gouraud", rasterized=True)
        ax.triplot(tri, color="0.2", lw=0.3, alpha=0.55, zorder=3)
        ax.axvline(0.0, color="0.4", ls="--", lw=1.3, zorder=5)
        F._style_axes(ax, None, ylabel=(key == "sv"), xlabel=True)
        fig.colorbar(tpc, ax=ax, location="bottom", fraction=0.09, pad=0.04, aspect=30).set_label(lab)

    axc = axd["cont"]
    contours(axc, tri, f, "solid", True)
    contours(axc, tri_m, pp, "dashed", False)
    F._style_axes(axc, None)
    stamp = f"openCSMP {yrs / 1000:g} kyrs, MPFA {F.snapshot_time_label(F.MPFA_DIR, None)}"
    axc.text(0.97, 0.95, stamp, transform=axc.transAxes, ha="right", va="top", fontsize=11)
    handles = [Line2D([], [], color=F.P_COLOR, lw=1.6, label="Pressure"),
               Line2D([], [], color="0.2", lw=1.6, label="openCSMP"),
               Line2D([], [], color=F.T_COLOR, lw=1.6, label="Temperature"),
               Line2D([], [], color="0.2", lw=1.6, ls="dashed", label="HU-PorePy MPFA"),
               Line2D([], [], color=F.S_COLOR, lw=1.6, label="Liquid saturation"),
               Line2D([], [], color="none", label=""),
               Line2D([], [], color=F.H_COLOR, lw=1.6, label="Enthalpy"),
               Line2D([], [], color="none", label="")]
    axc.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=4,
               fontsize=11, frameon=True, columnspacing=1.4, handlelength=1.6)

    out = args.out or f"csmp_snapshot_{yrs}yr"
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(F.HERE, f"{out}.{ext}"), dpi=300, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
