"""Copy the PorePy VTU snapshots behind the paper figures to benchmark_figures_data/porepy_vtks/.

Each entry below names a run folder written by the run scripts and the snapshot times the
paper uses (None = the final snapshot). The selected snapshots of every subdomain and mortar
grid are copied to porepy_vtks/<paper images folder>/<case>/, together with a data.pvd that
lists them with their times. Missing run folders are reported and skipped.

Usage:  python collect_porepy_vtks.py {4_1,4_2} [--list]
"""

import argparse
import os
import re
import shutil
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "benchmark_figures_data", "porepy_vtks")
DAY = 86400.0

_IM = "immiscible_multiphase_two_dimensional_flow"
_W1 = "weis_one_dimensional_flow"
_W2 = "weis_two_dimensional_flow"
_CO2 = "co2_two_dimensional_flow"
_1D = "subsection_4_2/one_dimensional"
_2D = "subsection_4_2/two_dimensional"

# (run folder, destination, snapshot times in days or None for the final one)
ENTRIES = {
    "4_1": [
        ("subsection_4_1/visualization_barriers_hu_N4", f"{_IM}/barriers_fixed_dim_N4", [0.0, 78.0, 571.0]),
        ("subsection_4_1/visualization_barriers_frac_hu_N3", f"{_IM}/barriers_mixed_dim_N3", [0.0, 78.0, 571.0]),
    ],
    "4_2": [
        *[(f"{_1D}/single_phase_visualization/case_{c}P_{o}", f"{_W1}/fig4_{c.lower()}p_{o}", None)
          for c in ("h", "m", "l") for o in ("horizontal", "vertical")],
        *[(f"{_1D}/visualization/fig5_{o}_hu", f"{_W1}/fig5_{o}_hu", None)
          for o in ("horizontal", "vertical")],
        *[(f"{_1D}/visualization/fig6_{c}", f"{_W1}/fig6_{c}", None) for c in ("pw", "salt")],
        (f"{_2D}/visualization_hu_simplex_trunc_q9", f"{_W2}/plume_tpfa", None),
        (f"{_2D}/visualization_hu_mpfa_simplex_trunc_q9", f"{_W2}/plume_mpfa", None),
        (f"{_2D}/visualization_co2_inject_hu_md", f"{_CO2}/co2_tpfa", None),
        (f"{_2D}/visualization_co2_inject_hu_md_mpfa", f"{_CO2}/co2_mpfa", None),
        (f"{_2D}/visualization_co2_inject_hu_md_mpfa_nograv", f"{_CO2}/co2_mpfa_no_gravity", None),
    ],
}


def _datasets(folder):
    """(time [s], snapshot index, vtu name) of every dataset in the folder's master pvd."""
    masters = [f for f in os.listdir(folder)
               if f.endswith(".pvd") and not re.search(r"_\d{6}\.pvd$", f)]
    if len(masters) != 1:
        raise FileNotFoundError(f"expected one master pvd in {folder}, found {masters}")
    out = []
    for ds in ET.parse(os.path.join(folder, masters[0])).getroot().iter("DataSet"):
        name = ds.attrib["file"]
        out.append((float(ds.attrib["timestep"]), int(re.search(r"_(\d{6})\.vtu$", name).group(1)), name))
    return out


def _select(datasets, days):
    """Snapshot indices nearest to the requested days (the last one for None)."""
    times = {idx: t for t, idx, _ in datasets}
    if days is None:
        return [max(times, key=times.get)]
    return [min(times, key=lambda i: abs(times[i] - d * DAY)) for d in days]


def collect(section, list_only=False):
    for src, dst, days in ENTRIES[section]:
        src_dir, dst_dir = os.path.join(HERE, src), os.path.join(OUT, dst)
        when = "final" if days is None else ", ".join(f"{d:g} d" for d in days)
        print(f"+ {src} ({when}) -> benchmark_figures_data/porepy_vtks/{dst}")
        if list_only:
            continue
        if not os.path.isdir(src_dir):
            print(f"  missing: {src}")
            continue
        datasets = _datasets(src_dir)
        keep = _select(datasets, days)
        shutil.rmtree(dst_dir, ignore_errors=True)
        os.makedirs(dst_dir)
        lines = []
        for t, idx, name in datasets:
            if idx in keep:
                shutil.copy2(os.path.join(src_dir, name), dst_dir)
                lines.append(f'\t<DataSet group="" part="" timestep="{t:.6f}" file="{name}"/>')
        with open(os.path.join(dst_dir, "data.pvd"), "w") as fh:
            fh.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" '
                     'byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n<Collection>\n'
                     + "\n".join(lines) + "\n</Collection>\n</VTKFile>\n")
        print(f"  copied {len(lines)} vtu files, snapshots {keep}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("section", choices=sorted(ENTRIES))
    ap.add_argument("--list", action="store_true", help="print what would be copied")
    args = ap.parse_args()
    collect(args.section, args.list)
