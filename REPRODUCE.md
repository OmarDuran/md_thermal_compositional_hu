# Reproducing the results

Run everything in Docker with a 16 GB memory limit, one script at a time.

## 1. Get the image

From the Zenodo archive:

```bash
docker load -i md_thermal_compositional_hu.tar.gz
```

Or build it from this repository, after placing the tables listed in the README:

```bash
docker build -t md_thermal_compositional_hu .
```

## 2. Run a subsection

The figures are written to `results/` on the host, in the layout of the paper's `images/` folder.

```bash
mkdir -p results
```

**Section 4.1**

```bash
docker run --rm --memory=16g --memory-swap=16g \
    -v "$PWD/results:/workdir/porepy/src/porepy/examples/geothermal_flow/paper_figures" \
    md_thermal_compositional_hu ./run_subsection_4_1.sh
```

**Sections 4.2 and 4.3 and Appendix B**

```bash
docker run --rm --memory=16g --memory-swap=16g \
    -v "$PWD/results:/workdir/porepy/src/porepy/examples/geothermal_flow/paper_figures" \
    md_thermal_compositional_hu ./run_subsection_4_2.sh
```

To run selected steps, name them and end with `collect`, which copies the figures to `results/`:

```bash
docker run --rm --memory=16g --memory-swap=16g \
    -v "$PWD/results:/workdir/porepy/src/porepy/examples/geothermal_flow/paper_figures" \
    md_thermal_compositional_hu ./run_subsection_4_2.sh weis1d collect
```

On Apple silicon the image runs under emulation. Single steps take hours and a full script can
take days.

## 3. Steps and figures

| Paper content | Script | Step | File in `results/` |
|---|---|---|---|
| Section 4.1, monotonicity and reduction table | `run_subsection_4_1.sh` | `checks` | `immiscible_multiphase_two_dimensional_flow/monotonicity.pdf`, `degenerate_equivalence.tex` |
| Section 4.1, column segregation | `run_subsection_4_1.sh` | `column` | `immiscible_multiphase_two_dimensional_flow/column_segregation_passive_phase.pdf` |
| Section 4.1, barriers | `run_subsection_4_1.sh` | `barriers` | `immiscible_multiphase_two_dimensional_flow/comparison_saturation_maps_{fixed_dim_4_phases,mixed_dim_3_phases}.pdf` |
| Section 4.2 1D, Figure 14 | `run_subsection_4_2.sh` | `weis1d` | `weis_one_dimensional_flow/reproduce_weis_fig_5.pdf` |
| Appendix B, Figure B.16 | `run_subsection_4_2.sh` | `weis1d` | `weis_one_dimensional_flow/reproduce_weis_fig_4.pdf` |
| Appendix B, Figure B.17 | `run_subsection_4_2.sh` | `weis1d` | `weis_one_dimensional_flow/reproduce_weis_fig_6.pdf` |
| Section 4.2 2D, plume | `run_subsection_4_2.sh` | `plume` | `weis_two_dimensional_flow/fig10_gravity_comparison.pdf` |
| Section 4.3, water–CO2 | `run_subsection_4_2.sh` | `co2` | `co2_two_dimensional_flow/co2_rho_h_contours.pdf`, `co2_breakthrough.pdf` |
| Appendix B, water–CO2 diagrams | `run_subsection_4_2.sh` | `phz` | `co2_two_dimensional_flow/co2_phz_diagrams.pdf` |

The OBL tables are included in the image.

## 4. Retrieve images from a container

Use this when the run was started without the `results/` mount. Give the container a name and
leave out `--rm`:

```bash
docker run --name hu_run --memory=16g --memory-swap=16g \
    md_thermal_compositional_hu ./run_subsection_4_2.sh weis1d collect
```

Copy the figures out. This also works while the container is running:

```bash
docker cp hu_run:/workdir/porepy/src/porepy/examples/geothermal_flow/paper_figures ./paper_figures
```

Copy the VTU snapshots behind the figures in the same way:

```bash
docker cp hu_run:/workdir/porepy/src/porepy/examples/geothermal_flow/benchmark_figures_data/porepy_vtks ./porepy_vtks
```

Remove the container when you are done:

```bash
docker rm hu_run
```
