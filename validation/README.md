# Validation suite

Experiments that validate VILMA2 beyond the unit tests of `make check`: the
community benchmarks, the lateral-viscosity benchmark, and the last deglaciation
run through VILMA2 and, behind the same driver, through VILMA1. One script,
`run.sh`, stages and submits the runs, reduces their output and draws the
figures and tables. The suite was moved here from the repository of the model
description paper (`alex-robinson/paper-fastearth3d-experiments`).

[`RESULTS.md`](RESULTS.md) records the last pass: model commit, jobs, outcome.

## Experiments

| Stage | What | Runs |
|---|---|---|
| `A` | Community benchmarks, each a pass/fail driver from `tests/`: Love numbers (A1), Spada (2011) disc load (A2), Martinec (2018) sea-level cases A, C2, D3, E2, F1 (A3), polar motion and its coupling into the sea-level equation (A4), Weerdesteijn (2023) low-viscosity zone at ℓmax 512 (A5), 3-D path against 1-D path and rotational invariance (A6), restart (A7) | 13 jobs, seconds to minutes each |
| `B` | Lateral-viscosity benchmark (Klemann et al. setup note), `tests/bench_visc3d.f90`: tests A-i, B-i, C-i with heaviside and ramp forcing at ℓmax 128, to 100 kyr in five chained 20-kyr segments. `LMAX=64` and `LMAX=256` with `CASES="A B"` give the truncation series | A-i, B-i: minutes. C-i: about 20 h per forcing at 64 threads |
| `ref` | LGM reference state: the 26 ka GLAC-1D load held for 30 kyr against RTopo-2 (3-D viscosity, ℓmax 170). Its end state is the `i_eq = 2` reference of every deglaciation run | 1 job |
| `D` | Last deglaciation, 26 ka to present, VILMA2 at ℓmax 170: 3-D Bagge et al. (2021) viscosity with the toroidal field off (`D-V2-3D.l170.fe`, V7's control) and the radial mean of that field (`D-V2-1Db`) | 2 jobs |
| `V1` | The same two runs through the VILMA1 backend (`solver = "v1"`), on VILMA2's layer table and the same viscosity values | 2 jobs |
| `T` | V7, the production 3-D run (`D-V2-3D.l170.tor`): `D`'s 3-D run with the toroidal field of the 3-D elements on (`l_toroidal`, the model default), otherwise identical | 1 job, about 1.5× the cost of the run without it |
| `E` | Sensitivity of the 3-D run: ℓmax 64, 128, 256 and `visc_3d_tol` | 7 jobs; scripts only, not part of the last pass |

Deglaciation setup (`experiments/nml/blockD.nml.tmpl`): PREM geometry in 11
layers, forward Euler at the Maxwell stability limit, rotation on, centre-of-mass
frame, ice density 931 kg m⁻³, viscosity clamped to 10^19.5–10^30 Pa s, the
220-element radial mesh, lossless output, all set explicitly whatever the model
defaults. The toroidal field is on in the production 3-D run (`T`) and off in its
control (`D`); the 1Db run has no laterally 3-D element for it to act on. The
analysis names these runs once, in `analysis/Published.jl` (`RUN_V2_3D`, ...). Both backends start from the same LGM reference file, so neither needs a
spin-up. `experiments/nml/blockDV1.nml.tmpl` is the VILMA1 counterpart.

## Running

```bash
cd validation
./run.sh build              # model, benchmark drivers, VILMA1-backend binary (login node)
./run.sh A B ref            # submit; wait for the queue
./run.sh submit prep        # one SLURM job: export the LGM reference, write the shared viscosity files
./run.sh D V1               # submit (from a login node: staging records the model commit with git)
./run.sh submit reduce analyse      # one SLURM job: runs/ -> data/published/ -> figures and tables
```

`./run.sh` without arguments lists the stages. `SUBMIT=0` stages without
submitting; `CASES`, `FORCINGS`, `LMAX` select block B cases; `ONLY="A1 A3"`
selects block A rows; `WALL` and `OMP` override wall time and threads of the
deglaciation runs.

Requirements:

- the model configured for the machine (`configme`), `runme` with
  `.runme/config.toml`, and for `V1` a VILMA1 installation (`include/`,
  `lib/vega_pism.a`) and its static input files;
- Julia ≥ 1.12 with the suite's own environment (`validation/Project.toml`,
  `Manifest.toml`: CairoMakie and NCDatasets) instantiated:
  `julia --project=. -e 'using Pkg; Pkg.instantiate()'` from `validation/`;
- `data/input/bagge2021.nc`: Bagge et al. (2021) log10 viscosity on lon, lat, r
  (512 × 256 × 164, quantised to 0.01 dex), the file this repository vendored
  as `input/bagge2021.nc` until commit 2676879. It derives from
  `isostasy_data/earth_structure/viscosity/bagge2021.nc`. Set `V_VISC_BAGGE`
  to use another copy.

## Machines

`machines/<name>.sh` holds what differs between clusters: account, partition,
cores per node, whether jobs take a whole node, threads and wall time of the
deglaciation runs, the VILMA1 installation and input directory, and the mapping
from wall time to a runme queue. The machine is chosen from the host name, or by
`VALIDATION_MACHINE=<name>`. `albedo` (AWI) and `levante` (DKRZ) exist; a new
machine is one new file.

## Layout

```
run.sh                 entry point, one stage per argument
Project.toml, Manifest.toml   the Julia environment of the analysis
machines/              per-cluster settings
experiments/           launch scripts, namelist templates (nml/), preparation of
                       the reference state and the shared viscosity files (*.jl),
                       vilma1_input_vilma2grid/densi.inp (VILMA2's layers for VILMA1)
analysis/              make_reduced.jl (the only reader of runs/), Published.jl
                       (the reader of data/published/), one script per block
tables/                generated tables and headline numbers (tracked)
RESULTS.md             record of the last pass
runs/                  raw output (untracked)
data/input/            bagge2021.nc (untracked)
data/reference/        LGM reference and derived viscosity files (untracked,
                       except visko_bagge.inp)
data/published/        reduced output the figures are made from (untracked)
figures/               diagnostic figures and all captions (untracked)
```

Figures shown in the documentation are written to `docs/benchmarks/figures/`
(tracked); `analysis/FigStyle.jl` lists them (`DOCS_FIGURES`).

On albedo the raw runs and the reduced data of the last pass are under
`/albedo/work/projects/p_forclima/robinson/models/vilma2-paper/validation/`
(`runs/`, `data/published/`, `data/reference/`).

## Notes

- `reduce` skips runs that are missing or incomplete, and `analyse` reads only
  `data/published/`. Re-run `reduce` after re-making a run.
- `D` and `V1` refuse to start until the reference run has written its own
  `spinup/vilma_restart.nc`; staging refuses a directory with a queued or
  running job.
- VILMA1 reads its viscosity file by level index and without a clamp. The files
  it gets are therefore written on its element centres with VILMA2's clamp
  applied (`make_vilma1_visc.jl`, `make_vilma1_visc_ongrid.jl`), and
  `check_vilma1_visc_grid.py` checks the levels before a run is staged.
- VILMA1 does not sub-step: `vilma_v1_nsub = 12` keeps its step (8.33 yr) below
  the shortest Maxwell time of the clamped field (10.5 yr).
- The VILMA1-backend binary is built with other objects than the default one;
  `build` therefore runs `make clean` twice. Do not run it in a tree another
  job is using.
