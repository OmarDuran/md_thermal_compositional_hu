# Weis et al. (2014) one-dimensional benchmarks

Figure 5 of Section 4.2 and Figures 4 and 6 of Appendix B. `../../run_subsection_4_2.sh weis1d`
runs this folder; the repository README lists all paper outputs.

| File | Purpose |
|---|---|
| `run_workflow.py` | one command for all three figures (resumable through `_cache/`) |
| `reproduce_weis_fig_{4,5,6}.py` | one figure each: single phase, two-phase boiling, brine with halite |
| `weis_1d_solver.py` | standalone 1-D finite-volume solver (PPU, HU, HU-mw), imported by the figure scripts |
| `porepy_1d_solver.py`, `single_phase_porepy_1d_solver.py` | PorePy runs for the HU overlays |
| `fig_weis_common.py`, `plot_style.py` | shared caching and plot style |
| `fig6_{pw,salt}_ph_consistency.py` | Fig. 6 diagnostics, not in the paper |

```bash
python run_workflow.py --serial --pdf   # full run, writes figures/reproduce_weis_fig_{4,5,6}.{png,pdf}
python run_workflow.py --plot-only      # figures from the existing _cache/ only
python run_workflow.py --no-porepy      # standalone solver only, no PorePy overlays
python run_workflow.py --quick          # short test run in _quick/, real cache untouched
```

Delete files in `_cache/` to recompute those runs.
