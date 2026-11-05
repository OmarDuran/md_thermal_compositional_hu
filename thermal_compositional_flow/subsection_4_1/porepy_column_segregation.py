#!/usr/bin/env python
"""Column segregation with a passive intermediate phase (paper Figure 7), PorePy only.

A vertical column of height 100 m, density-inverted: heavy alpha on top, light gamma at the
bottom, a passive intermediate phase beta in the central third. Three densities of beta sweep
the alpha-gamma pair,

    near-heavy   rho_beta = 1250   (m_beta = 1)
    exact mean   rho_beta = 1000   (m_beta = 1/2)
    near-light   rho_beta =  750   (m_beta = 0)

and each is run with HU (``porepy_2d_solver --scheme hu``) and PPU (``--scheme ppu``) on the
same model, grid, time stepping and PETSc linear solver. The model is the subsection_4_1
barriers solver with the geometry and initial bands overridden.

Run:  python porepy_column_segregation.py                 # dy = 1 m, 60 d, figure
      python porepy_column_segregation.py --ny 200        # dy = 0.5 m
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import porepy as pp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import porepy_2d_solver as base  # noqa: E402

RHO_A, RHO_G = 1500.0, 500.0  # heavy (alpha) / light (gamma) [kg/m^3]
CASES = {  # label -> rho_beta
    "near_heavy": 0.75 * RHO_A + 0.25 * RHO_G,  # 1250
    "exact": 0.50 * RHO_A + 0.50 * RHO_G,  # 1000
    "near_light": 0.25 * RHO_A + 0.75 * RHO_G,  # 750
}
CASE_TITLE = {
    "near_heavy": r"near-heavy, $\rho_\beta = 1250$",
    "exact": r"exact mean, $\rho_\beta = 1000$",
    "near_light": r"near-light, $\rho_\beta = 750$",
}
SCHEMES = ("hu", "ppu")
DISPLAY = {"hu": "HU", "ppu": "PPU"}

LY = 100.0  # column height [m]
WIDTH = 1.0  # column width [m] (nx = 1)


def midpoint_weight(rho_beta: float) -> float:
    """Density weight m_beta of the (alpha, gamma) pair, in {0, 1/2, 1}."""
    s = np.sign((RHO_A - RHO_G) * (2.0 * rho_beta - RHO_A - RHO_G))
    return 0.5 * (1.0 + s)


class _ColumnOverrides(pp.PorePyModel):
    """Barrier-free 1 x ny column with three equal bands; records s_beta(y, t) per step."""

    _height = LY
    _width = WIDTH

    def set_domain(self) -> None:
        H = self.units.convert_units(self._height, "m")
        W = self.units.convert_units(self._width, "m")
        self._domain = pp.Domain({"xmax": W, "ymax": H})

    def grid_type(self):
        return "cartesian"

    def meshing_arguments(self) -> dict:
        H = self.units.convert_units(self._height, "m")
        W = self.units.convert_units(self._width, "m")
        ny = int(self.params.get("ny", 100))
        return {"cell_size": H / ny, "cell_size_x": W}

    def set_fractures(self) -> None:
        self._fractures = []

    def barrier_cell_mask(self, sd: pp.Grid) -> np.ndarray:
        return np.zeros(sd.num_cells, dtype=bool)

    def _band_bounds(self) -> list:
        H = self.units.convert_units(self._height, "m")
        return list(np.linspace(H, 0.0, base.NPHASE + 1))

    def initial_condition(self) -> None:
        super().initial_condition()
        self._seg_t, self._seg_prof, self._seg_y = [], [], None
        self._record_profile()

    def after_nonlinear_convergence(self) -> None:
        super().after_nonlinear_convergence()
        self._record_profile()

    def _record_profile(self) -> None:
        sd = self.mdg.subdomains(dim=self.nd)[0]
        order = np.argsort(sd.cell_centers[1])  # bottom -> top
        if self._seg_y is None:
            self._seg_y = sd.cell_centers[1][order]
        beta = list(self.fluid.phases)[1]  # phase 1 = intermediate (beta)
        s_beta = np.asarray(beta.saturation([sd]).value(self.equation_system), dtype=float)
        self._seg_t.append(float(self.time_manager.time) / base.day)
        self._seg_prof.append(s_beta[order])


class ColumnFlowModel(_ColumnOverrides, base.FlowModel):
    """HU on the column."""


class ColumnPPUFlowModel(_ColumnOverrides, base.PPUFlowModel):
    """PPU on the column."""


def run_suffix(constant_dt: bool, dt_days: float, nl_tol: float) -> str:
    """File-name tag for non-default runs (empty for the default adaptive, tol 1e-4 run)."""
    tag = f"_cdt{dt_days:g}" if constant_dt else ""
    return tag + (f"_tol{nl_tol:.0e}" if nl_tol != 1.0e-4 else "")


def run_case(rho_beta: float, scheme: str, *, ny: int = 100, perm_factor: float = 10.0,
             t_end_days: float = 60.0, dt_init_days: float = 0.05, dt_max_days: float = 0.25,
             out_root: str = "visualization_column", constant_dt: bool = False,
             nl_tol: float = 1.0e-4) -> dict:
    """Run one column to ``t_end_days``; return the s_beta(y, t) history and solver counts."""
    k_rock = base.k_rock * perm_factor
    solid = pp.SolidConstants(
        permeability=k_rock, porosity=base.porosity, thermal_conductivity=2.0 * base.to_Mega,
        density=2500.0, specific_heat_capacity=base.C_P,
        residual_aperture=base.FRACTURE_APERTURE,
    )
    tag = f"{scheme}_{int(round(rho_beta))}_ny{ny}" + run_suffix(constant_dt, dt_max_days, nl_tol)
    params = base.build_params(
        3, scheme, t_end_days=t_end_days, dt_init_days=dt_init_days, dt_max_days=dt_max_days,
        snap_days=(0.0, t_end_days), fractures=False, ny=ny,
        material_constants={"solid": solid},
        folder_name=os.path.join(out_root, f"col_{tag}"), file_name=f"col_{tag}",
    )
    if constant_dt:
        # Ramp from dt_init (doubling) to dt_max, then hold dt_max: dt only shrinks on a failed step.
        day = 86400.0
        params["time_manager"] = pp.TimeManager(
            schedule=[0.0, t_end_days * day], dt_init=dt_init_days * day, constant_dt=False,
            dt_min_max=(dt_max_days * day / 64.0, dt_max_days * day),
            iter_optimal_range=(13, 13), iter_relax_factors=(0.5, 2.0),
            recomp_factor=0.5, recomp_max=8)
    # The middle-phase density of this case (read by the EOS and the z -> s inversion).
    base.RHO = np.array([RHO_A, float(rho_beta), RHO_G])
    base.RHO_BANDS = base.RHO.copy()

    model_class = ColumnPPUFlowModel if scheme == "ppu" else ColumnFlowModel
    model = model_class(params)
    n_steps = max(1, round(t_end_days / dt_max_days))
    drift_tol = 10.0 ** (-(base.DRIFT_ORDER - 1)) / (2 * n_steps)
    solver_params = {
        "nl_convergence_criteria": {
            "res_abs": pp.solvers.ResidualBasedAbsoluteCriterion(
                tol=nl_tol, metric=pp.EquationBasedLebesgueMetric(model)),
            "null_drift": base.NullSpaceDriftCriterion(model, tol=drift_tol),
        },
        "nl_divergence_criteria": {
            "max_iter": pp.solvers.MaxIterationsCriterion(max_iterations=13),
        },
    }
    runner = pp.ModelRunner(model, solver_params,
                            nonlinear_solver=base.geothermal_nonlinear_solver(solver_params))
    runner.run()

    stats = model.collect_run_stats()
    t = np.asarray(model._seg_t)
    y = np.asarray(model._seg_y)
    prof = np.asarray(model._seg_prof)  # (ntime, ny), bottom -> top
    w = np.clip(prof, 0.0, None)
    com = (w * y[None, :]).sum(axis=1) / np.maximum(w.sum(axis=1), 1e-30)
    out = dict(t=t, y=y, prof=prof, com=com, m_beta=midpoint_weight(rho_beta),
               iters=int(stats.total_newton_iterations), steps=int(stats.n_accepted_steps),
               cuts=int(stats.n_time_step_cuts))
    os.makedirs(out_root, exist_ok=True)
    np.savez(os.path.join(out_root, f"{tag}.npz"), **out)
    print(f"  {scheme:3s} rho_b={rho_beta:6.0f}  com {com[0]:.2f} -> {com[-1]:.2f} m  "
          f"steps={out['steps']} newton={out['iters']} cuts={out['cuts']}", flush=True)
    return out


def make_figure(results: dict, out_dir: str, ny: int, t_end_days: float,
                perm_factor: float) -> None:
    """Top: HU s_beta(y, t) with the HU and PPU centre-of-mass tracks.
    Bottom: |s_beta^HU - s_beta^PPU|, PPU interpolated to the HU times."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        import seaborn as sns

        cmap = sns.color_palette("vlag", as_cmap=True)
    except Exception:  # pragma: no cover
        cmap = plt.get_cmap("coolwarm")

    labels = list(CASES)
    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.6), constrained_layout=True,
                             sharex=True, sharey=True)
    im_s = im_d = None
    for k, lab in enumerate(labels):
        hu, ppu = results[lab]["hu"], results[lab]["ppu"]
        ax, axd = axes[0, k], axes[1, k]
        im_s = ax.pcolormesh(hu["t"], hu["y"], hu["prof"].T, cmap=cmap, vmin=0.0, vmax=1.0,
                             shading="gouraud", rasterized=True)
        ax.plot(hu["t"], hu["com"], color="k", lw=1.6,
                label=f"HU ({hu['iters']} Newton)")
        ax.plot(ppu["t"], ppu["com"], color="#15b01a", lw=1.5, ls="--",
                label=f"PPU ({ppu['iters']} Newton)")
        ax.set_title(f"{CASE_TITLE[lab]}\n$m_\\beta = {hu['m_beta']:.1f}$", fontsize=12)
        ax.legend(loc="lower right", fontsize=9, framealpha=0.85)

        # PPU profiles on the HU time levels (the two adaptive step sequences differ).
        ppu_on_hu = np.array([[np.interp(tt, ppu["t"], ppu["prof"][:, j]) for j in
                               range(ppu["prof"].shape[1])] for tt in hu["t"]])
        diff = np.abs(hu["prof"] - ppu_on_hu)
        im_d = axd.pcolormesh(hu["t"], hu["y"], diff.T, cmap="magma_r", vmin=0.0,
                              vmax=max(0.05, float(diff.max())), shading="gouraud",
                              rasterized=True)
        axd.set_title(f"max $|s_\\beta^{{HU}} - s_\\beta^{{PPU}}| = {diff.max():.3f}$",
                      fontsize=11)
        axd.set_xlabel("time [d]")
        if k == 0:
            ax.set_ylabel("height $y$ [m]")
            axd.set_ylabel("height $y$ [m]")
        ax.set_xlim(0, t_end_days)
        ax.set_ylim(0, LY)

    fig.colorbar(im_s, ax=list(axes[0]), location="right", fraction=0.04, pad=0.02,
                 shrink=0.9).set_label(r"$s_\beta$ (HU)")
    fig.colorbar(im_d, ax=list(axes[1]), location="right", fraction=0.04, pad=0.02,
                 shrink=0.9).set_label(r"$|s_\beta^{HU} - s_\beta^{PPU}|$")
    fig.suptitle(
        f"PorePy HU and PPU column segregation, $\\Delta y = {LY / ny:g}$ m, "
        f"$K = {perm_factor * 1000.0 / 1000.0:g}$ D, {t_end_days:.0f} d", fontsize=13)
    os.makedirs(out_dir, exist_ok=True)
    stem = f"porepy_column_segregation_hu_ppu_ny{ny}"
    for ext in ("png", "pdf"):
        path = os.path.join(out_dir, f"{stem}.{ext}")
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print("wrote", os.path.relpath(path, HERE))
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ny", type=int, default=100, help="column cells (default 100: dy = 1 m)")
    ap.add_argument("--perm-factor", type=float, default=10.0,
                    help="rock permeability factor (default 10: K = 10 D)")
    ap.add_argument("--t-end-days", type=float, default=60.0)
    ap.add_argument("--dt-init-days", type=float, default=0.05)
    ap.add_argument("--dt-max-days", type=float, default=0.25)
    ap.add_argument("--cases", nargs="+", default=list(CASES), choices=list(CASES))
    ap.add_argument("--schemes", nargs="+", default=list(SCHEMES), choices=list(SCHEMES))
    ap.add_argument("--out", default=os.path.join(HERE, "passive_phase_example"))
    ap.add_argument("--constant-dt", action="store_true", help="ramp from --dt-init-days to --dt-max-days, then hold it")
    ap.add_argument("--tol", type=float, default=1.0e-4, help="Newton residual tolerance")
    args = ap.parse_args()

    print(f"=== PorePy column segregation: ny={args.ny} (dy={LY / args.ny:g} m), "
          f"K x{args.perm_factor:g}, {args.t_end_days:g} d ===", flush=True)
    results: dict = {}
    for lab in args.cases:
        results[lab] = {}
        for scheme in args.schemes:
            results[lab][scheme] = run_case(
                CASES[lab], scheme, ny=args.ny, perm_factor=args.perm_factor,
                t_end_days=args.t_end_days, dt_init_days=args.dt_init_days,
                dt_max_days=args.dt_max_days, constant_dt=args.constant_dt, nl_tol=args.tol)
    if (set(args.schemes) == set(SCHEMES) and set(args.cases) == set(CASES)
            and not run_suffix(args.constant_dt, args.dt_max_days, args.tol)):
        make_figure(results, args.out, args.ny, args.t_end_days, args.perm_factor)


if __name__ == "__main__":
    main()
