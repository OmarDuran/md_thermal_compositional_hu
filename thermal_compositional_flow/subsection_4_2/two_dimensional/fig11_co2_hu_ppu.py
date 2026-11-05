#!/usr/bin/env python
"""Water-CO2 figure with HU and PPU (gravity, MPFA).

(a) rows HU and |HU - PPU|; columns mixture density and mixture enthalpy, last snapshot.
(b) outlet CO2 fraction with HU, one curve per discretization.

Reads the run folders (read-only). Writes fig11_co2_hu_ppu_{maps,breakthrough}.{png,pdf} to --out.

Run:  python fig11_co2_hu_ppu.py --out DIR
"""
import argparse
import os

import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402

import plot_rho_h_contours as C  # noqa: E402  (colormap, ranges, snapshot lookup)

HERE = os.path.dirname(os.path.abspath(__file__))
HU, PPU = "visualization_co2_inject_hu_md_mpfa", "visualization_co2_inject_ppu_md_mpfa"
RUNS = [  # (label, folder, colour, linestyle); HU only
    ("TPFA", "visualization_co2_inject_hu_md", "C0", "-"),
    ("MPFA", "visualization_co2_inject_hu_md_mpfa", "C1", "-"),
    ("MPFA, no gravity", "visualization_co2_inject_hu_md_mpfa_nograv", "C2", "-"),
]


def _fields(folder, field):
    f2, _ = C._latest_file(os.path.join(HERE, folder), 2)
    m2 = pv.read(f2)
    v2 = np.asarray(m2.cell_data[field], float)
    f1, _ = C._latest_file(os.path.join(HERE, folder), 1)
    m1 = pv.read(f1) if f1 else None
    v1 = np.asarray(m1.cell_data[field], float) if m1 is not None and field in m1.cell_data else None
    return m2, v2, m1, v1


def _draw(ax, m2, v2, m1, v1, vmin, vmax, cmap):
    cc = m2.cell_centers().points
    ax.tricontourf(cc[:, 0], cc[:, 1], v2, levels=np.linspace(vmin, vmax, C.NLEV), cmap=cmap,
                   vmin=vmin, vmax=vmax, extend="both")
    e = m2.extract_all_edges()
    conn = e.lines.reshape(-1, 3)[:, 1:]
    ax.add_collection(LineCollection(e.points[conn][:, :, :2], colors="0.35", linewidths=0.15,
                                     zorder=2))
    if m1 is not None and v1 is not None:
        segs = m1.points[m1.cells_dict[3]][:, :, :2]
        ax.add_collection(LineCollection(segs, colors="0.15", linewidths=C.FRAC_LW + 1.2,
                                         capstyle="projecting", zorder=4))
        lc = LineCollection(segs, cmap=cmap, norm=plt.Normalize(vmin, vmax),
                            linewidths=C.FRAC_LW, capstyle="projecting", zorder=5)
        lc.set_array(v1)
        ax.add_collection(lc)
    ax.set_xlim(0, C.L_DOMAIN)
    ax.set_ylim(0, C.L_DOMAIN)
    ax.set_aspect("equal")


def maps(out):
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 11.5), sharex=True, sharey=True,
                             gridspec_kw=dict(wspace=0.20, hspace=0.06))
    for c, (field, title, clabel, (vmin, vmax), cmap) in enumerate(C.COLS):
        hu = _fields(HU, field)
        ppu = _fields(PPU, field)
        _draw(axes[0, c], *hu, vmin, vmax, cmap)
        d2 = np.abs(hu[1] - ppu[1])
        d1 = np.abs(hu[3] - ppu[3]) if hu[3] is not None and ppu[3] is not None else None
        dmax = float(d2.max()) if d2.max() > 0 else 0.01 * (vmax - vmin)
        _draw(axes[1, c], hu[0], d2, hu[2], d1, 0.0, dmax, cmap)
        axes[0, c].set_title(title)
        axes[1, c].set_xlabel("x [m]")
        sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin, vmax), cmap=cmap)
        fig.colorbar(sm, ax=axes[0, c], location="bottom", shrink=0.9, pad=0.04).set_label(clabel)
        sm = plt.cm.ScalarMappable(norm=plt.Normalize(0.0, dmax), cmap=cmap)
        fig.colorbar(sm, ax=axes[1, c], location="bottom", shrink=0.9, pad=0.16).set_label(
            "|HU $-$ PPU| " + clabel.split(" ", 1)[1])
    for r, lab in enumerate(("HU", "difference")):
        axes[r, 0].set_ylabel(f"{lab}\ny [m]")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"fig11_co2_hu_ppu_maps.{ext}"), dpi=150, bbox_inches="tight")


def breakthrough(out):
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    for label, folder, col, ls in RUNS:
        path = os.path.join(HERE, folder, "outlet_co2_breakthrough.csv")
        if not os.path.exists(path):
            continue
        d = np.genfromtxt(path, delimiter=",", names=True)
        t, z = np.atleast_1d(d["time_days"]), np.atleast_1d(d["outlet_z_co2_avg"])
        ax.plot(t, z, color=col, ls=ls, lw=1.6, label=label)
        k = np.flatnonzero(z > 1e-3)
        print(f"{label:18s} first z>1e-3 at {t[k[0]] if k.size else float('nan'):6.2f} d, "
              f"z(end={t[-1]:.1f} d) = {z[-1]:.4f}")
    ax.set_xlabel("time [days]")
    ax.set_ylabel(r"outlet CO$_2$ fraction [-]")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"fig11_co2_hu_ppu_breakthrough.{ext}"), dpi=150,
                    bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=HERE)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    maps(a.out)
    breakthrough(a.out)


if __name__ == "__main__":
    main()
