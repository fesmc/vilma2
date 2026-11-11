# VILMA2

**VILMA2** is a state-of-the-art but simple and fast **3D solid-Earth model** — a
visco-elastic deformation model coupled with the sea-level equation — intended
as an **open-source alternative to VILMA1** within the CLIMBER-X climate model.

The method is a new implementation of the **spectral–finite-element,
time-domain** approach of Martinec (2000): spherical harmonics horizontally,
finite elements radially, an incompressible Maxwell rheology integrated
explicitly in time, a self-consistent sea-level equation with migrating
coastlines, and rotational feedback. It is built **3D-ready from the start**
(laterally varying viscosity) and validated against the published GIA benchmarks.

**Documentation:** <https://fesmc.github.io/vilma2/> (physics,
discretization, benchmarks, install & run), rendered from the Quarto sources
under [`docs/`](docs/); see [doc/design.md](doc/design.md) for the design
rationale and method comparison, and
[doc/vilma-v1-backend.md](doc/vilma-v1-backend.md) for the optional VILMA1 backend
(`&vilma solver = "v1"`), which drives VILMA1 itself through this model's
driver, namelist, forcing, remap and output for a like-for-like comparison. It
is **off by default and is not a dependency**: it needs an explicit
`make vilma vilma_v1=1 VILMA_V1_ROOT=<install>` and a hand-installed VILMA1.

## Status

VILMA2 is validated against the published GIA benchmarks and is starting to be
used for science. It is versioned 0.x: the namelist and the host API may still
change. Each benchmark is an executable test with a pass/fail exit
(`make check`, `make check-slow`); the reference build is ifx `-Ofast` on the
AWI cluster albedo.

| Benchmark | Test | Asserted | Measured (albedo, ifx) |
|---|---|---|---|
| Love numbers, M3-L70-V01 (Spada et al. 2011) | `test_benchmark_love` | 1 %: elastic h, k at degrees 2–48; fluid at 2–8 | h ≤ 0.08 % to degree 256; k, l exceed 1 % above degrees 96 and 144 |
| Disc load (Spada et al. 2011) | `test_benchmark_disc` | elastic profile 6 % / 4 % of peak; center uplift 3 % at 1–10 kyr | 4.0 % / 2.8 %; ≤ 0.6 % |
| Cap load, case A (Martinec et al. 2018) | `test_benchmark_martinec` | 1.5 % (u, horizontal), 3 % (geoid) | 0.5 %, 0.7 %, 1.9 % |
| Sea-level equation, cases C2, D3, E2, F1 (Martinec et al. 2018) | `test_benchmark_sle` | 30 % (u), 35 % (horizontal), 12 % (geoid, rsl), 8 % (sea surface) of peak | D3, E2, F1 ≤ 7.6 %, rsl ≤ 1.8 %; C2 up to 28 % (unattributed degree-0/1 offset) |
| Polar motion, Spada test 3/2 | `test_rotation` | 3 % | ≤ 0.9 % |
| Low-viscosity zone (Weerdesteijn et al. 2023) | `test_benchmark_lvz` | 8 % | 2.6 %, 1.0 % |

The model builds and passes the test suite with ifx on albedo and with gfortran
(macOS, Ubuntu CI); it builds with ifx on DKRZ levante. Open items: the C2
sea-level case, the high-degree k and l Love numbers (limited by the radial
mesh, second order in its spacing), and the 3-D viscosity
path, which is compared only with VILMA1 and the axisymmetric
low-viscosity-zone case. The toroidal flow driven by lateral viscosity is
implemented but **off by default** (`l_toroidal = .false.`) for speed: it makes
3-D runs 1.6–1.8 times slower. A production deglaciation run with it on is
planned to quantify its effect. Details:
[benchmarks](https://fesmc.github.io/vilma2/benchmarks/).

## Install

Dependencies are [fesm-utils](https://github.com/fesmc/fesm-utils) at branch
`phase0-kryos` (FFTW, SHTns, the `fesmutils` helper library and the `coords`
module) plus a system netCDF (C and Fortran). The branch is pinned in
`.configme/manifest.toml`: VILMA2 needs `coords`' `grid_init_gaussian`, which
fesm-utils `dev` (v1.4) does not have yet. **configme** clones and builds
fesm-utils and generates the machine/compiler Makefile:

```bash
configme install vilma2 -m macbook -c gfortran      # or -m awi_albedo / dkrz_levante -c ifx
configme install vilma2 --link fesm-utils=/abs/path/to/fesm-utils   # reuse a checkout
```

The exact sequence for the macbook, albedo and levante is in
[Install & run](https://fesmc.github.io/vilma2/getting-started.html).

## Build

```bash
make vilma        # -> bin/vilma.x        (standalone forced-run driver)
make vilma_mkref  # -> bin/vilma_mkref.x  (build a Gauss-grid reference)
make vilma_remap  # -> bin/vilma_remap.x  (offline lon-lat -> Gauss remap)
make vilma-static # -> libvilma/include/libvilma.a and the .mod files
make check                 # build + run the test suite (about 2 min)
make openmp=1 check-slow   # full-resolution benchmarks and order tests (about 10 min on 32 threads)
```

Both test targets run every test, end with a summary naming the failures and
exit non-zero if any failed. Setting `VILMA_BENCH_DUMP=<dir>` makes the
benchmark tests write their comparison tables there.

`make` switches: `debug=0|1|2`, `openmp=0|1` (default 1; build `openmp=1` for the
threaded degree loop at production resolutions).

### As a library

`make vilma-static` writes the objects, the `.mod` files and `libvilma.a` to
`libvilma/include/` (as yelmo and chion). The archive is a file target: its
time changes only when an object does, so a host can depend on it. `openmp=0`
and `openmp=1` share the directory; a stamp (`.openmp0|1`) rebuilds all objects
when the setting changes, and the host must use the same setting as the
fesm-utils variant it links (`include-serial` or `include-omp`). A host
compiles with `-I<vilma>/libvilma/include` and links, in this order,

```
-L<vilma>/libvilma/include -lvilma  -L<fesm-utils>/include-omp -lfesmutils
-L<fesm-utils>/SHTns/shtns-omp/lib -lshtns_omp
-L<fesm-utils>/fftw/fftw-omp/lib -lfftw3_omp -lfftw3 -lm   <netCDF>
```

(`-serial`, `-lshtns`, `-lfftw3` for `openmp=0`).

### Host-resolved terms (optional)

A host that resolves the ice margin on its own grid inside regions of weight `w` can pass, per
unit Gauss-cell area, the grounded-ice column `I_g = w·H·f_g`, the ocean fraction
`C = w·(1−f_g)` and the flooded-area fill `F = w·ζ⁽⁰⁾·(C − C_ref)`:

```fortran
call solid_earth_init(se, z_bed_eq, h_ice_eq, grid, w=w, I_g=I_g, C=C)    ! reference terms
call solid_earth_spinup(se, h_ice_lgm, w, I_g, C, F)                      ! held over every pass
call solid_earth_update(se, h_ice, dt_yr, w, I_g, C, F)                   ! ramped like the ice
ok = solid_earth_host_check(w, I_g, C, F, msg)                            ! admissibility
```

Inside `w` the sea-level solve uses them; on `1 − w` it applies its own flotation and coastline
rule, and `h_ice`, `h_ice_eq` are then the ice of that part (an area mean; cells with `w = 1`
ignore them). All or none, fixed for the run; not with `solver = "v1"` or `sle_fixed_ocean`;
`F` enters only through the subgrid term (`sle_subgrid`). Without them, or with
`w ≡ 0`, a run is bit-identical to one without. Fields on the Gauss grid pass through; a global
lon-lat host (`grid`) is remapped by the model, a regional one is refused. The terms are part of
the restart. Details: [doc/design.md §15](doc/design.md).

## Configure & run

A run configuration has two namelist groups: `&vilma` (physics and numerics,
the `vilma_param_class` record a host fills in memory) and `&ctl` (standalone
run control: forcing, reference and output files, time window, `i_eq`,
restart). The complete, documented `&vilma` defaults are in
[`input/vilma_defaults.nml`](input/vilma_defaults.nml) and are loaded
automatically; a run config such as [`vilma.nml`](vilma.nml) (the default) or
[`examples/deglac_lgm.nml`](examples/deglac_lgm.nml) carries the `&ctl` group
plus only the `&vilma` parameters it overrides (yelmo `defaults_file`
convention). Time fields (`dt_*`, `time_*`) are given in **years** and
converted to SI seconds on load.

- **Earth structure** — `earth`: a named built-in (`"M3-L70-V01"`, `"PREM"`) or
  `"custom"` to assemble from the surface-first layer arrays.
- **Response solver** — `earth_response`: `"ve"` (full viscoelastic, default),
  `"elastic"`, `"null"`.
- **Memory scheme** — `scheme = "fe"` (1st-order explicit, default) or `"trap"`
  (2nd-order adaptive), advanced by the `vilma_timestep` controller.
- **Degree-1 frame** — `deg1_frame = "cm"` (default, geocenter motion kept) or
  `"cf"` (Spada disc-benchmark convention).
- **3D viscosity** — `l_visc_3d`, `visc_3d_file`, `visc3d_tol`; `l_toroidal`
  carries the toroidal flow (off by default). `visc_3d_structure = n` reads
  structure `n` of a FesmData-layout file instead of a lon-lat-r field. The
  default is the FesmData product
  `fesmdata/Global/GLOBAL-T170/GLOBAL-T170_SolidEarth-Bagge2021.nc` (Bagge et
  al. 2021, structure 6 = `v_1.0_s16`): the path is relative to the run folder,
  which needs a `fesmdata` link to the FesmData store (Levante:
  `/work/ba1442/fesmdata`). No viscosity file is shipped with the model.
- **Spin-up / restart** — `equil_time_max`, `pre_spinup_1d`, and
  `restart_in_file` (`&ctl`).
- **Rotation** — `rotation` (TPW feedback): on by default; `.false.` for the
  non-rotating benchmarks.
- **Backend** — `solver = "v2"` (default) or `"v1"` (see above).

Run the standalone driver on a run config:

```bash
./bin/vilma.x examples/deglac_lgm.nml    # default: vilma.nml
```

It reads a reference state and an ice-thickness forcing (`file_forcing`,
`h_ice(lon,lat,time)`; remapped from lon-lat on the fly by default), marches the
model across the forcing, and writes the diagnostic surface fields (`rsl`,
`z_bed`, …) to `file_out`, and optionally the horizontal displacement to
`file_hor`.

Or stage/submit runs and ensembles with **runme** (`-r` run, `-s` submit;
comma-lists in `-p` define ensemble dimensions):

```bash
runme -o runs/deglac -e main --omp 8 -r -p vilma.lmax=128 vilma.earth_response=ve
```

Embedding the model in a host (the CLIMBER-X coupling path) uses the same API
behind a single `use vilma`:

```fortran
use vilma
type(solid_earth) :: se
call vilma_par_load(se%par, "vilma.nml")
call solid_earth_init(se, z_bed_eq, h_ice_eq, grid=host_grid, time=t0)  ! Gauss grid + remap; clock [yr]
call solid_earth_update(se, h_ice, dt_yr)   ! advance dt_yr [years]; read se%rsl, se%z_bed
                                            ! (rsl = se%dz_ss - se%dz_bed: sea-surface and bed change)
call solid_earth_finalize(se)
```

`use vilma` is the only module a host uses (the contract, no internals). `h_ice` is the whole
ice column; the model applies its own flotation. The model owns its clock: `solid_earth_init`
sets it (`time`, default 0), updates advance it, `solid_earth_set_time` relabels it; no caller
writes `se%time`. A standalone run resumes a restart at the restart's time (`&ctl
restart_time_reset` overrides), and a restart from another mesh or viscosity is refused
(`&ctl restart_other_visc` accepts it as a seed).

## License

VILMA2 is free software under the GNU General Public License v3.0; see
[LICENSE](LICENSE).
