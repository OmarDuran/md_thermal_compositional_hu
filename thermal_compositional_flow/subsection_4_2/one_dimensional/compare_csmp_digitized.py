"""Compare the openCSMP 1D profiles (csmp_1d_solver.py) with the digitized Weis (2014) curves.

For every digitized curve the openCSMP profile is interpolated to the digitized points; the
table reports the RMS and largest differences, and one figure per case overlays both.

    python compare_csmp_digitized.py [cases...]
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE.parent.parent / "benchmark_figures_data"
RUNS = HERE / "csmp_visualization"

# digitized-file suffix -> (profile column, scale to the digitized unit, label, unit)
FIELDS = {
    "pressure": ("pressure_Pa", 1.0e-6, "pressure", "MPa"),
    "pressured": ("pressure_Pa", 1.0e-6, "pressure", "MPa"),
    "temperature": ("temperature_C", 1.0, "temperature", "degC"),
    "saturation_liq": ("saturation_liquid", 1.0, "liquid saturation", "-"),
    "saturation_halite": ("saturation_halite", 1.0, "halite saturation", "-"),
}


def read_digitized(path):
    with open(path) as f:
        rows = [r for r in csv.reader(f)][1:]
    xy = np.array([[float(a), float(b)] for r in rows if len(r) >= 2 and r[0].strip() and r[1].strip() for a, b in [r[:2]]])
    return xy[np.argsort(xy[:, 0])]


def read_profile(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def compare(case):
    prof_path = RUNS / case / "csmp_profile_1.csv"
    if not prof_path.exists():
        print(f"{case}: no openCSMP profile")
        return
    prof = read_profile(prof_path)
    x_km = prof["x_m"] / 1000.0
    base = case.split("_n")[0]                           # fig_5_vertical_n800 -> fig_5_vertical
    files = sorted(DATA.glob(f"{base}_*_raw.csv"))
    fig, axes = plt.subplots(1, len(files), figsize=(4.2 * len(files), 3.4))
    for ax, f in zip(np.atleast_1d(axes), files):
        field = f.name[len(base) + 1:-len("_raw.csv")]
        col, scale, label, unit = FIELDS[field]
        dig = read_digitized(f)
        sim = np.interp(dig[:, 0], x_km, prof[col] * scale)
        err = sim - dig[:, 1]
        span = np.ptp(dig[:, 1]) or 1.0
        print(f"{case:22s} {label:18s} RMS {np.sqrt(np.mean(err**2)):9.4g} {unit:4s} "
              f"max {np.max(np.abs(err)):9.4g} {unit:4s} ({100 * np.max(np.abs(err)) / span:5.1f} % of range)")
        ax.plot(x_km, prof[col] * scale, "-", color="C0", label="openCSMP")
        ax.plot(dig[:, 0], dig[:, 1], "o", mfc="none", color="k", ms=4, label="Weis (2014), digitized")
        ax.set_xlabel("distance [km]")
        ax.set_ylabel(f"{label} [{unit}]")
    np.atleast_1d(axes)[0].legend(fontsize=8)
    fig.suptitle(case)
    fig.tight_layout()
    out = RUNS / case / f"{case}_vs_digitized.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    cases = sys.argv[1:] or sorted(p.name for p in RUNS.iterdir() if p.is_dir())
    for c in cases:
        compare(c)
