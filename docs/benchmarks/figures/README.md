# Benchmark figures

The benchmark pages reference these files (made by the scripts in the
repository's `validation/` directory from the `VILMA_BENCH_DUMP` tables):

| File | Page |
|---|---|
| `love-numbers.png` | love-numbers.qmd |
| `disc.png` | disc.qmd |
| `martinec-A.png`, `sle.png`, `sle-C2.png`, `sle-D3.png`, `sle-E2.png`, `sle-F1.png` | sea-level.qmd |
| `rotation.png` | rotation.qmd |
| `lvz.png` | lvz.qmd |
| `visc3d-sections-heav.png`, `visc3d-sections-ramp.png`, `visc3d-geocentre.png`, `visc3d-A-normal-mode.png`, `visc3d-lmax.png` | visc3d.qmd (hidden page) |
| `deglac-maps.png`, `deglac-v2v1-maps-3d.png`, `deglac-v2v1-maps-1d.png`, `deglac-v2v1-summary.png` | deglaciation.qmd (hidden page) |

`validation/RESULTS.md` records the model commit and the jobs of the pass that
made them. To regenerate them: `cd validation && ./run.sh submit reduce analyse`
after the runs (see `validation/README.md`).

A page that references a figure not yet made can hold its image line in a
comment block, `<!-- FIGURE (enable when figures/NAME.png exists) ... -->`, so
that the site renders without it; the line below enables all such blocks, from
the `docs/` directory:

```bash
perl -0pi -e 's/<!-- FIGURE \(enable when [^)]*\)\n(.*?)\n-->/$1/gs' benchmarks/*.qmd
```
