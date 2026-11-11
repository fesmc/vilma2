# Benchmark figures

The benchmark pages reference these files (made by the scripts in the
repository's `validation/` directory from the `VILMA_BENCH_DUMP` tables):

| File | Page |
|---|---|
| `love-numbers.png` | love-numbers.qmd |
| `disc.png` | disc.qmd |
| `sle.png`, `sle-C2.png`, `sle-D3.png`, `sle-E2.png`, `sle-F1.png` | sea-level.qmd |
| `rotation.png` | rotation.qmd |
| `lvz.png` | lvz.qmd |
| `visc3d-sections-heav.png`, `visc3d-sections-ramp.png`, `visc3d-geocentre.png`, `visc3d-A-normal-mode.png`, `visc3d-lmax.png` | visc3d.qmd (hidden page) |
| `deglac-v2v1-maps-3d.png`, `deglac-v2v1-maps-1d.png`, `deglac-v2v1-summary.png` | deglaciation.qmd (hidden page) |

Until a file is committed here, its image line sits in a comment block of the
form

```
<!-- FIGURE (enable when figures/NAME.png exists)
![caption](figures/NAME.png){#fig-NAME}
-->
```

so that the site renders without it. Once the files are here, enable all blocks
in one edit from the `docs/` directory:

```bash
perl -0pi -e 's/<!-- FIGURE \(enable when [^)]*\)\n(.*?)\n-->/$1/gs' benchmarks/*.qmd
```
