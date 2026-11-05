# md_thermal_compositional_hu

Scripts reproducing the numerical results of "A Mixed-Dimensional Compositional Flow Model
Incorporating Gravitational Effects". The model is non-isothermal compositional flow of $N$
phases in fractured porous media, written in fractional flow form, with buoyancy discretized by
hybrid upwinding (HU) and phase properties evaluated by operator-based linearization (OBL). The
scripts extend the `geothermal_flow` example of [PorePy](https://github.com/pmgbergen/porepy),
branch `buoyancy_flow_refactor`.

## Contents

| Step | Paper content | Output in `images/` |
|---|---|---|
| `run_subsection_4_1.sh checks` | monotonicity of the buoyant flux; reduction consistency for $N = 3, 4$ | `immiscible_multiphase_two_dimensional_flow/monotonicity.pdf`, `degenerate_equivalence.tex` |
| `run_subsection_4_1.sh column` | passive-phase column, $\rho_\beta \in \{750, 1000, 1250\}$ kg/m³ | `.../column_segregation_passive_phase.pdf` |
| `run_subsection_4_1.sh barriers` | $N$-phase segregation through barriers, fixed- and mixed-dimensional, HU and PPU | `.../comparison_saturation_maps_*.pdf` |
| `run_subsection_4_2.sh weis1d` | H₂O–NaCl columns of Weis et al. (2014), Figures 10, B.13 | `weis_one_dimensional_flow/reproduce_weis_fig_{5,4,6}.pdf` |
| `run_subsection_4_2.sh plume` | 2D plume, TPFA versus MPFA gravity, Figure 11 | `weis_two_dimensional_flow/fig10_gravity_comparison.pdf` |
| `run_subsection_4_2.sh co2` | three-phase H₂O–CO₂ across conductive and blocking fractures, Figure 12 | `co2_two_dimensional_flow/co2_{rho_h_contours,breakthrough}.pdf` |
| `run_subsection_4_2.sh phz` | H₂O–CO₂ phase diagrams, Figure B.15 | `co2_two_dimensional_flow/co2_phz_diagrams.pdf` |

The final `collect` step copies the figures to `thermal_compositional_flow/paper_figures/` and the
PorePy snapshots behind them to `benchmark_figures_data/porepy_vtks/`. The openCSMP references
come from `run_opencsmp.sh` (1D columns and the plume of Figure 11; `save` packs all its output
into `opencsmp_data.tar.gz`), which must run before `run_subsection_4_2.sh weis1d`. The H₂O–NaCl
error maps of Appendix B are produced by the table builder in `opensowat`.

## Running

The image is archived on Zenodo, https://zenodo.org/records/23872069. Load it, or build it from
this repository.

```bash
docker load -i md_thermal_compositional_hu.tar.gz    # or: docker build -t md_thermal_compositional_hu .
docker run --rm --memory=16g --memory-swap=16g \
    -v "$PWD/results:/workdir/porepy/src/porepy/examples/geothermal_flow/paper_figures" \
    md_thermal_compositional_hu ./run_subsection_4_2.sh weis1d collect
```

Run one script at a time within 16 GB; `--list` prints the commands of a script. The image is
linux/amd64 and runs under emulation on Apple silicon. `REPRODUCE.md` gives the full steps and the
retrieval of results from a container.

The build needs files too large for git on disk: the OBL tables listed in `.dockerignore` and the
openCSMP lookup tables in `opencsmp/tables/`.

## License

GPL-3.0. openCSMP (LGPL-3.0) is built in the image from a pinned commit and
`opencsmp/opencsmp_weis.patch`, with its source kept in `/workdir/openCSMP`.
