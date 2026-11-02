"""openCSMP (CSMP++) runs of the one-dimensional Weis et al. (2014) benchmarks, Figs. 4-6.

Each case is a 2 km column with fixed pressure, temperature and salinity at both ends and a
linear initial pressure (condition 1). The openCSMP driver ``csmp_1d_column`` (built in the
openCSMP repository) runs it with the CVFEM PHX scheme on 200 elements (dx = 10 m, as in the
paper); this script writes the configuration and run files, runs the cases one after another
and leaves the final profile as ``csmp_visualization/<case>/csmp_profile_1.csv``.

    python csmp_1d_solver.py                      # all ten cases
    python csmp_1d_solver.py fig_5_vertical       # selected cases
    python csmp_1d_solver.py --list
    python csmp_1d_solver.py --elements 800 fig_5_vertical   # PorePy resolution, own folder
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "csmp_visualization"
YEAR = 365.0 * 86400.0

LENGTH = 2000.0           # [m]
ELEMENTS = 200            # dx = 10 m (Weis et al. 2014)
DT_INITIAL = 0.01         # [yr]
DT_MAX = 1.0              # [yr], the paper's largest step
DT_MIN = 1.0e-8           # [yr], failure guard
CFL_SCALING = 0.1
ROCK = dict(permeability=1.0e-15, porosity=0.1, thermal_conductivity=2.0,
            density=2700.0, specific_heat_capacity=880.0, compressibility=1.0e-15)

# name: (p_left, p_right [MPa], T_left, T_right [degC], gravity, final time [yr], extras);
# x = 0 is the inflow end and, with gravity, the bottom
CASES = {
    "fig_4_hp_horizontal": (50.0, 25.0, 350.0, 150.0, False, 250.0, {}),
    "fig_4_hp_vertical":   (50.0, 25.0, 350.0, 150.0, True, 750.0, {}),
    "fig_4_mp_horizontal": (40.0, 20.0, 450.0, 300.0, False, 120.0, {}),
    "fig_4_mp_vertical":   (40.0, 20.0, 450.0, 300.0, True, 350.0, {}),
    "fig_4_lp_horizontal": (15.0, 1.0, 500.0, 350.0, False, 1500.0, {}),
    "fig_4_lp_vertical":   (15.0, 1.0, 500.0, 350.0, True, 1500.0, {}),
    "fig_5_horizontal":    (20.0, 1.0, 400.0, 150.0, False, 200.0, {}),
    "fig_5_vertical":      (20.0, 1.0, 400.0, 150.0, True, 1000.0, {}),
    "fig_6_pw":            (4.0, 1.0, 300.0, 150.0, False, 2000.0, {}),
    # salt-saturated liquid with halite saturation 0.1, pure water at the inflow end, and Weis
    # eq. 28 for the halite permeability (Weis et al. 2014, Fig. 6C,D)
    "fig_6_salt":          (4.0, 1.0, 300.0, 150.0, False, 2000.0,
                            dict(target_halite=0.1, halite_permeability=1, salinity_right=40.0)),
}


def write_config(path: Path, case: tuple) -> None:
    """openCSMP configuration: blocks are separated by exactly one blank line."""
    p_l, p_r, t_l, t_r, gravity, tf, extra = case
    p_l, p_r = p_l * 1.0e6, p_r * 1.0e6
    w_r = extra.get("salinity_right", 0.0)
    cap = ROCK["specific_heat_capacity"] * ROCK["density"] * (1.0 - ROCK["porosity"])
    k = ROCK["permeability"]
    path.write_text(f"""# (block1) region name from parameter range

# (block2) default property values
fluid pressure\t\t\t{p_r:g}
temperature\t\t\t{t_r:g}
permeability\t\t\t{k:g}
vertical permeability\t\t{k:g}
horizontal permeability\t\t{k:g}
porosity\t\t\t{ROCK['porosity']:g}
nodal porosity\t\t\t{ROCK['porosity']:g}
salinity\t\t\t{w_r:g}
thermal conductivity\t\t{ROCK['thermal_conductivity']:g}
nodal heat capacity rock\t{ROCK['specific_heat_capacity']:g}
total heat capacity\t\t{cap:g}
nodal heat capacity\t\t{cap:g}
density rock\t\t\t{ROCK['density']:g}
nodal density rock\t\t{ROCK['density']:g}
lithostatic pressure\t\t{p_r:g}
nodal compressibility rock\t{ROCK['compressibility']:g}
gravity flag\t\t\t{1.0 if gravity else 0.0:g}
cfl scaling\t\t\t{CFL_SCALING:g}
open top\t\t\t0
no standard output\t\t1
thickness\t\t\t1.

# (block3) regional property values

# (block4) boundary conditions for box-shaped model
LEFT\tDIRICH\tfluid pressure\t{p_l:g}\t{p_l:g}
RIGHT\tDIRICH\tfluid pressure\t{p_r:g}\t{p_r:g}
LEFT\tDIRICH\ttemperature\t{t_l:g}\t{t_l:g}
RIGHT\tDIRICH\ttemperature\t{t_r:g}\t{t_r:g}
LEFT\tDIRICH\tsalinity\t0.\t0.
RIGHT\tDIRICH\tsalinity\t{w_r:g}\t{w_r:g}

# (block5) regional property conditions

# (block6) boundary conditions for arbitrary-shaped model

# (block7) run settings [s]
duration\t{tf * YEAR:g}
time increment\t{DT_MAX * YEAR:g}
output time\t{tf * YEAR:g}

""")


def find_driver(explicit: str | None) -> Path:
    candidates = [explicit, os.environ.get("CSMP_1D_COLUMN")]
    candidates += [str(p / "openCSMP" / "build" / "csmp_1d_column") for p in HERE.parents]
    for c in candidates:
        if c and Path(c).is_file() and os.access(c, os.X_OK):
            return Path(c).resolve()
    raise SystemExit("csmp_1d_column not found: build it (ninja -C build csmp_1d_column in openCSMP) "
                     "and pass --csmp-exe or set $CSMP_1D_COLUMN")


def link_tables(tables: Path, run_dir: Path) -> int:
    n = 0
    for tab in tables.glob("*LookupTable*.bin"):
        target = run_dir / tab.name
        if not target.exists():
            target.symlink_to(tab.resolve())
        n += 1
    return n


def run_case(name: str, driver: Path, tables: Path, elements: int = ELEMENTS) -> int:
    case = CASES[name]
    run_dir = OUT / (name if elements == ELEMENTS else f"{name}_n{elements}")
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "run_complete.json").unlink(missing_ok=True)
    write_config(run_dir / "csmp-configuration.txt", case)
    extra = case[6]
    run = [f"config {run_dir / 'csmp'}", f"output {run_dir / 'csmp'}",
           f"length {LENGTH:g}", f"elements {elements}", f"snap_years 0 {case[5]:g}",
           f"dt_initial {DT_INITIAL:g}", f"dt_max {DT_MAX:g}", f"dt_min {DT_MIN:g}",
           "linear_rtol 1e-10", "linear_maxit 1000", f"cfl_scaling {CFL_SCALING:g}",
           f"rock_permeability {ROCK['permeability']:g}"]
    if "target_halite" in extra:
        run.append(f"target_halite {extra['target_halite']:g}")
    if extra.get("halite_permeability"):
        run.append("halite_permeability 1")
    run_file = run_dir / "csmp.run"
    run_file.write_text("\n".join(run) + "\n")
    if link_tables(tables, run_dir) == 0:
        print(f"NOTE: no lookup tables in {tables}; openCSMP builds them in this run (about an hour).")

    print(f"=== {name}: {case[0]:g} -> {case[1]:g} MPa, {case[2]:g} -> {case[3]:g} degC, "
          f"{'vertical' if case[4] else 'horizontal'}, {case[5]:g} yr, {elements} elements", flush=True)
    t0 = time.time()
    with open(run_dir / "csmp.log", "w") as out:
        proc = subprocess.Popen([str(driver), str(run_file)], cwd=run_dir, stdout=out,
                                stderr=subprocess.PIPE, text=True, bufsize=1)
        markers = ("csmp_1d_column", "column:", "initial pressure", "initial halite")
        for line in proc.stderr:
            out.write(line)
            hits = [line.find(m) for m in markers if m in line]
            if hits:
                print(line[min(hits):].rstrip(), flush=True)
        status = proc.wait()
    elapsed = time.time() - t0
    if status == 0:
        with open(run_dir / "run_complete.json", "w") as f:
            json.dump({"completed": True, "final_time_years": case[5], "elapsed_run_seconds": elapsed,
                       "elements": elements}, f, indent=2)
        print(f"{name}: done in {elapsed:.0f} s", flush=True)
    else:
        print(f"{name}: FAILED with status {status}; see {run_dir / 'csmp.log'}", flush=True)
    return status


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cases", nargs="*", help="case names (default: all)")
    ap.add_argument("--list", action="store_true", help="list the cases and stop")
    ap.add_argument("--csmp-exe", help="path to csmp_1d_column")
    ap.add_argument("--elements", type=int, default=ELEMENTS,
                    help=f"elements over the 2 km column (default {ELEMENTS}; PorePy uses 800)")
    ap.add_argument("--tables-dir", help="openCSMP lookup-table cache (default: next to the driver)")
    args = ap.parse_args()
    if args.list:
        print("\n".join(CASES))
        return 0
    names = args.cases or list(CASES)
    unknown = [n for n in names if n not in CASES]
    if unknown:
        raise SystemExit(f"unknown cases {unknown}; use --list")
    driver = find_driver(args.csmp_exe)
    tables = Path(args.tables_dir or os.environ.get("CSMP_TABLES_DIR")
                  or driver.parent / "example_inputs" / "CVFEM_examples")
    print(f"openCSMP driver: {driver}\noutput folder:   {OUT}")
    failed = [n for n in names if run_case(n, driver, tables, args.elements) != 0]
    if failed:
        print(f"failed: {failed}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
