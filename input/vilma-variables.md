# vilma restart / output variables

Variable-io table (yelmo convention) for VILMA2 netCDF I/O, used by
`vilma_io` for both restart files and diagnostic `write_step` output. Each row gives
the netCDF variable name, its dimensions, units, and a long_name. The time axis
(unlimited) lets several snapshots live in one file.

The prognostic state restored on restart depends on the response kind: the
Maxwell memory-stress fields `tau_*` (+ `sigma_n_*`) for `RESP_VE`, plus the
adaptive controller's next-step Δt
seed `dt_try_s` (`dt_try` is the same in years, for reading) for both. The reference fields `z_bed_eq`/`h_ice_eq` are static
(written once) and checked on read, and so are the host reference terms
`w_ref`, `I_g_ref`, `C_host_ref` of a run with host-resolved terms; the host terms
of the last update (`w_prev`, `I_g_prev`, `C_prev`, `F_prev`: the next ramp start)
are prognostic; the rest are diagnostic.

The `prec` column is the output precision of the variable (ncio bit rounding, parsed by
`variable_io`): `-` lossless, `1e-3` absolute (units of the field), `1e-4r` relative,
`1e-3,1e-4r` both. It applies to diagnostic output only: restart files are always written
lossless, and `write_lossless` (`&ctl`, or `nc_set_lossless`) forces lossless output too.
Diagnostic surface fields are rounded to 1e-3 m (horizontal displacements 1e-4 m); the
prognostic state, the masks and the time series stay lossless.

| id | variable     | dimensions    | units  | long_name                                        | prec      |
|----|--------------|---------------|--------|--------------------------------------------------|-----------|
|  1 | tau_a_re     | nlam, ne, nk  | Pa     | Maxwell memory stress, component A (real part)   | -         |
|  2 | tau_a_im     | nlam, ne, nk  | Pa     | Maxwell memory stress, component A (imag part)   | -         |
|  3 | tau_b_re     | nlam, ne, nk  | Pa     | Maxwell memory stress, component B (real part)   | -         |
|  4 | tau_b_im     | nlam, ne, nk  | Pa     | Maxwell memory stress, component B (imag part)   | -         |
|  5 | tau_c_re     | nlam, ne, nk  | Pa     | Maxwell memory stress, component C (real part)   | -         |
|  6 | tau_c_im     | nlam, ne, nk  | Pa     | Maxwell memory stress, component C (imag part)   | -         |
|  7 | z_bed_eq     | lon, lat      | m      | Reference (equilibrium) bedrock elevation        | -         |
|  8 | h_ice_eq     | lon, lat      | m      | Reference (equilibrium) grounded-ice thickness   | -         |
|  9 | h_ice        | lon, lat      | m      | Grounded-ice thickness                           | 1e-3      |
| 10 | rsl          | lon, lat      | m      | Relative sea-level change (full field)           | 1e-3      |
| 11 | z_bed        | lon, lat      | m      | Bedrock elevation (z_bed_eq - rsl)               | 1e-3      |
| 12 | C_ocean      | lon, lat      | 1      | Ocean function (1 ocean / 0 land)                | -         |
| 13 | dt_try       | time          | years  | Adaptive time-stepping next-step Δt suggestion   | -         |
| 14 | sigma_n_re   | nlm, time     | kg m-2 | Trapezoidal start-of-step load σ_n (real part)   | -         |
| 15 | sigma_n_im   | nlm, time     | kg m-2 | Trapezoidal start-of-step load σ_n (imag part)   | -         |
| 16 | sigma_primed | time          | 1      | Flag: σ_n is tracked (1) or not yet primed (0)   | -         |
| 17 | bsl          | time          | m      | Barystatic sea level vs reference (ocean mean of C*rsl) | -         |
| 18 | dz_bed       | lon, lat      | m      | Bed change u (solid uplift); rsl = dz_ss - dz_bed | 1e-3      |
| 19 | dz_ss        | lon, lat      | m      | Sea-surface change N + dphi (geoid rise + shift) | 1e-3      |
| 20 | rot_m_re     | time          | 1      | Polar motion m₁ (rotation, real part)            | -         |
| 21 | rot_m_im     | time          | 1      | Polar motion m₂ (rotation, imag part)            | -         |
| 22 | rot_load_mem | nlam_rot, ne_rot, nrc | Pa | Rotation loading-channel memory stress (packed) | -         |
| 23 | rot_tidal_mem| nlam_rot, ne_rot, nrc | Pa | Rotation tidal-channel memory stress (packed)   | -         |
| 24 | u_east       | lon, lat      | m      | Surface horizontal displacement, east (total)    | 1e-4      |
| 25 | u_north      | lon, lat      | m      | Surface horizontal displacement, north (total)   | 1e-4      |
| 26 | u_east_tor   | lon, lat      | m      | Surface horizontal displacement, east (toroidal) | 1e-4      |
| 27 | u_north_tor  | lon, lat      | m      | Surface horizontal displacement, north (toroidal)| 1e-4      |
| 28 | dphi         | time          | m      | Uniform sea-surface shift (part of dz_ss)        | -         |
| 29 | bsl_ive      | time          | m      | Ice-volume equivalent of the grounded-ice change | -         |
| 30 | w_ref        | lon, lat      | 1      | Host reference: fraction of the cell in host regions | -         |
| 31 | I_g_ref      | lon, lat      | m      | Host reference: w H f_g, grounded-ice column     | -         |
| 32 | C_host_ref   | lon, lat      | 1      | Host reference: w (1 - f_g), ocean fraction      | -         |
| 33 | w_prev       | lon, lat      | 1      | Host terms of the last update: w (ramp start)    | -         |
| 34 | I_g_prev     | lon, lat      | m      | Host terms of the last update: I_g (ramp start)  | -         |
| 35 | C_prev       | lon, lat      | 1      | Host terms of the last update: C (ramp start)    | -         |
| 36 | F_prev       | lon, lat      | m      | Host terms of the last update: F (ramp start)    | -         |
| 37 | fill_int     | time          | m3     | Flooded-area fill of the last update, whole sphere: the term of the ocean mass balance | -         |
| 38 | fill_int_own | time          | m3     | Flooded-area fill, own part's share (a zeta0 (C_own - C0_own)) | -         |
| 39 | dt_try_s     | time          | s      | Adaptive time-stepping next-step Δt suggestion, in seconds (exact; dt_try is the same in years) | -         |
