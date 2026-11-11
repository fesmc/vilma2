"""
    Published

The reader for `data/published/` — the **reduced, publishable** model output
that every figure in this repository is drawn from.

WHY THIS EXISTS. The raw output is tens of GB of `runs/` on a cluster, untracked.
A reader has none of it, and neither does the author on a laptop. So
the pipeline is cut in two:

    runs/  ──(analysis/make_reduced.jl)──▶  data/published/  ──(plot scripts)──▶  figures/
           only this step touches runs/                     only this module reads it

`data/published/` is a few tens of megabytes of NetCDF and travels with the
paper. Once it exists, every figure can be redrawn — and the analysis behind it
modified — with no cluster, no 850 MB `out.nc`, and no model.

THE CONTRACT, which is what makes the split safe:

  * **The plot scripts never open `runs/`.** They call this module and nothing
    else. If `data/published/` is missing they fail immediately with a message
    naming the stage that builds it, rather than half-working off whatever raw
    output happens to be lying around.
  * **This module never opens `runs/`.** It cannot fall back, by construction.
  * **The schema is declared here** — `EPOCHS`, `SITES`, `V2V1_SERIES_TIMES`, the
    file names — and `make_reduced.jl` includes this file to get them. One
    declaration, obeyed by the writer and the reader, so the two cannot drift.
  * **Anything that cannot be rebuilt from the reduced fields is precomputed**
    during reduction and stored as a series. That is everything whose derivation
    sweeps all 261 time levels: the sea-level budget, the V2−V1 residual series,
    the raw-field series, the per-site time series.

WHAT IS AND IS NOT AVAILABLE OFFLINE. Maps are stored at the `EPOCHS` below and
nowhere else, so `slice_at(r, "rsl", -16000)` raises rather than silently
returning the 14 ka field. Asking for an epoch that was not reduced is a request
to re-run `make_reduced.jl` with a longer `EPOCHS`, on the machine that still
has `runs/`. Everything on the full `time` axis — barystatic sea level, polar
motion, the site series, the budget, the V2−V1 residual — is complete.

ORIENTATION. Fields are stored the way they are plotted, not the way the model
wrote them: longitude rolled to −180…180 (so the Atlantic is not split down the
middle of the frame) and latitude ascending south→north. The roll-and-flip that
`BlockD.orient` applies on every read is therefore done once, at reduction, and
the published files are directly usable by anything that reads NetCDF.
"""
module Published

using NCDatasets
using Printf
using Statistics

export published_dir, have_published, require_published,
       EPOCHS, SITES, V2V1_SERIES_TIMES, BUDGET_ROWS,
       PubRun, load_run, load_runs_blockD,
       slice_at, rsl_anomaly, uplift_rate, srf_elevation, nearest_epoch, nearest_t,
       ocean_mean, polar_motion, site_series, regrid, lat_weights,
       budget_series, budget_at, raw_series, deg01,
       load_v2v1, load_blockE, ka, sh_parts,
       BenchRun, load_blockA, dumptable, hasdump, verdict, passed, linear_nodes,
       BENCH_TESTS, BENCH_FORCINGS, BENCH_EPOCHS_KYR, BENCH_NPT, BENCH_LMAX, BENCH_CONV_LMAX,
       BENCH_CONV_TESTS, BenchV3D, load_blockB, bench_love_path, bench_cost

const REPO = dirname(@__DIR__)

published_dir() = joinpath(REPO, "data", "published")
blockD_dir()    = joinpath(published_dir(), "blockD")
v2v1_path()       = joinpath(published_dir(), "v2v1.nc")
blockE_path()   = joinpath(published_dir(), "blockE.nc")
blockA_path()   = joinpath(published_dir(), "blockA.nc")
blockB_dir()    = joinpath(published_dir(), "bench_visc3d")

# ---------------------------------------------------------------------------
# The schema
# ---------------------------------------------------------------------------
"""
The epochs at which two-dimensional fields are kept, in model years (negative =
past). Everything the figures draw as a map is drawn at one of these.

  −21000  the LGM panel, and the epoch every V2−V1 and 3-D−1Db table is quoted at
  −14000  the mid-deglaciation panel of the lateral-signal figure
  −10000  the mid-deglaciation panel of the deglaciation-maps figure
   −8000  the second row of the V2−V1 map figures
   −2000  the last row of the V2−V1 tables
    −500  ONLY to difference against 0 for the present-day uplift rate; it is
          not plotted anywhere on its own
       0  present day — the datum every anomaly is measured against, so it is
          required by every map whether or not it is drawn

Adding an epoch costs about 2.8 MB per run before compression and requires a
re-run of `make_reduced.jl` on the machine that holds `runs/`.
"""
const EPOCHS = [-21000.0, -14000.0, -10000.0, -8000.0, -2000.0, -500.0, 0.0]

"""
How far an epoch request may miss a stored epoch before it is an error, in
years. One coupling interval. The stored epochs are exact time levels of the
model's own 100 yr axis, so a legitimate request lands on one exactly; anything
further away is asking for a field that was not reduced, and must say so rather
than quietly return the neighbour.
"""
const EPOCH_TOL = 100.0

"""
The sites whose full relative-sea-level history is carried. Near field (Hudson
Bay, Ångermanälven), intermediate (Iceland, Patagonia), far field (Barbados,
Tahiti), Antarctic (Amundsen, Ross). `(name, lat, lon)` with lon in −180…180.

This list is part of the published schema rather than of the figure that plots
it: the series are extracted during reduction, so a site added here is not
available until `make_reduced.jl` runs again.
"""
const SITES = [
    ("Hudson Bay",     60.0,  -85.0),
    ("Ångermanälven",  63.0,   17.8),
    ("Iceland",        64.9,  -19.0),
    ("Barbados",       13.1,  -59.6),
    ("Tahiti",        -17.6, -149.6),
    ("Patagonia",     -50.0,  -73.0),
    ("Amundsen Sea",  -75.0, -105.0),
    ("Ross Sea",      -78.0,  175.0),
]

"""
The epochs at which the VILMA2−VILMA1 residual is decomposed. Every entry
needs both runs' full fields differenced and a degree-0/1 fit, so this is
computed once at reduction and stored as a series.
"""
const V2V1_SERIES_TIMES = collect(-26000.0:1000.0:0.0)

"Number of time levels the sea-level budget is evaluated at, spread over the run."
const BUDGET_ROWS = 27

have_published() = isdir(published_dir())

"""
    require_published(what = "")

Fail with the command that builds the published data, instead of with a
`No such file` from three frames down. Called by every loader here.
"""
function require_published(path::AbstractString)
    isfile(path) || isdir(path) || error("""
        missing published data: $(relpath(path, REPO))

        The figures are drawn from data/published/, not from runs/. Build it with

            ./run.sh reduce

        on a machine that holds runs/ (this repository's raw output is not
        distributed), or obtain data/published/ alongside the paper and unpack it
        here.
        """)
    path
end

# ---------------------------------------------------------------------------
# Block D — the deglaciation runs
# ---------------------------------------------------------------------------
"""
    PubRun

One reduced deglaciation run. `lon`/`lat` are the oriented axes, `epoch` the
stored map epochs and `time` the full coupling axis the series live on.

`ds` is held open: the field arrays are read on demand, which keeps a script
that draws one epoch from paging in all seven.
"""
struct PubRun
    id::String
    path::String
    lmax::Int
    visc::String          # "1d" | "3d"
    solver::String        # "VILMA2" | "VILMA1"
    lon::Vector{Float64}  # −180…180, ascending
    lat::Vector{Float64}  # south→north, ascending
    epoch::Vector{Float64}
    time::Vector{Float64}
    ds::NCDataset
end

Base.show(io::IO, r::PubRun) =
    print(io, "PubRun($(r.id), ℓmax $(r.lmax), $(r.visc), $(length(r.time)) steps, ",
              "$(length(r.epoch)) epochs)")

"""
    load_run(id) -> PubRun | nothing

`load_run("D-V2-3D.l170.fe")`. `nothing` when that run was not reduced, so a
caller can skip a case the way it used to skip a missing run directory.
"""
function load_run(id::AbstractString)
    require_published(blockD_dir())
    path = joinpath(blockD_dir(), String(id) * ".nc")
    isfile(path) || return nothing
    ds = NCDataset(path)
    PubRun(String(id), path,
           Int(ds.attrib["lmax"]), String(ds.attrib["visc"]), String(ds.attrib["solver"]),
           Float64.(ds["lon"][:]), Float64.(ds["lat"][:]),
           Float64.(ds["epoch"][:]), Float64.(ds["time"][:]), ds)
end

"Every reduced block D run, by id."
function load_runs_blockD()
    require_published(blockD_dir())
    d = Dict{String,PubRun}()
    for f in sort(readdir(blockD_dir()))
        endswith(f, ".nc") || continue
        id = f[1:end-3]
        r = load_run(id)
        r === nothing || (d[id] = r)
    end
    d
end

"""
    nearest_epoch(r, yr) -> Int

Index into `r.epoch`. Raises when `yr` is further than `EPOCH_TOL` from every
stored epoch — the reduced files carry seven map epochs and no others, and
returning the nearest one regardless would answer a question that was not asked.
"""
function nearest_epoch(r::PubRun, yr::Real)
    i = argmin(abs.(r.epoch .- yr))
    if abs(r.epoch[i] - yr) > EPOCH_TOL
        error("""
            $(r.id): no reduced field at $(yr) yr.
            Stored epochs: $(join(string.(Int.(r.epoch)), ", ")) yr.
            Add it to Published.EPOCHS and re-run analysis/make_reduced.jl on a
            machine that holds runs/.
            """)
    end
    i
end

"""
    nearest_t(r, yr) -> Int

Index into `r.time`, the FULL coupling axis the published series live on. This
is not `nearest_epoch`, and the difference matters: the series — barystatic sea
level, polar motion, the site histories — are complete at every one of the 261
levels, while the two-dimensional fields exist only at the seven `EPOCHS`. A
caller wanting `bsl` at 21 ka uses this; a caller wanting a map at 21 ka uses
`nearest_epoch`, which will refuse an epoch that was not reduced.
"""
nearest_t(r::PubRun, yr::Real) = argmin(abs.(r.time .- yr))

"""
    slice_at(r, name, yr) -> Matrix (lon, lat)

One stored epoch of a two-dimensional field. `z_bed` is derived rather than
stored: the model defines it as `z_bed_eq - rsl` exactly, so keeping it would be
a quarter of the file for no information.
"""
function slice_at(r::PubRun, name::AbstractString, yr::Real)
    name == "z_bed" && return Float64.(r.ds["z_bed_eq"][:, :]) .- slice_at(r, "rsl", yr)
    Float64.(r.ds[name][:, :, nearest_epoch(r, yr)])
end

"""
    rsl_anomaly(r, yr) -> Matrix [m]

Relative sea level at `yr` minus the model's own present-day field. Positive =
sea level higher then than now.
"""
rsl_anomaly(r::PubRun, yr::Real) = slice_at(r, "rsl", yr) .- slice_at(r, "rsl", 0.0)

"""
    uplift_rate(r; yr = 0.0, window = 500.0) -> Matrix [mm/yr]

Bedrock uplift rate, differenced over `window` years. `z_bed = z_bed_eq − rsl`,
so this is −Δrsl/Δt, positive up. The 500 yr window is why `−500` is in `EPOCHS`.
"""
function uplift_rate(r::PubRun; yr::Real = 0.0, window::Real = 500.0)
    i1 = nearest_epoch(r, yr)
    i0 = nearest_epoch(r, yr - window)
    dt = r.epoch[i1] - r.epoch[i0]
    dt == 0 && error("uplift_rate: zero time window at yr=$yr")
    a1 = Float64.(r.ds["rsl"][:, :, i1])
    a0 = Float64.(r.ds["rsl"][:, :, i0])
    @. -(a1 - a0) / dt * 1000.0
end

"""
    srf_elevation(r, yr; h_min = 10.0) -> (z_srf, ice_only)

Surface elevation `z_bed + h_ice`, and the same masked to NaN where the ice is
thinner than `h_min`. NaN renders transparent in Makie, which is how the ice
overlay sits on top of the sea-level field without hiding it.
"""
function srf_elevation(r::PubRun, yr::Real; h_min::Real = 10.0)
    h = slice_at(r, "h_ice", yr)
    z = slice_at(r, "z_bed", yr)
    zs = z .+ h
    ice = copy(zs)
    ice[h .< h_min] .= NaN
    (zs, ice)
end

"""
    lat_weights(lat) -> Vector

Quadrature weights in sin(lat), normalised to 1, for either latitude ordering.
The model grid is Gaussian, so plain cos(lat) weights are close but not exact —
close enough to hide a real error in a budget quoted to a tenth of a metre.
"""
function lat_weights(lat)
    n = length(lat); s = sind.(Float64.(lat)); e = zeros(n + 1)
    descending = s[1] > s[end]
    e[1]     = descending ?  1.0 : -1.0
    e[n + 1] = descending ? -1.0 :  1.0
    for j in 2:n
        e[j] = 0.5 * (s[j - 1] + s[j])
    end
    w = [abs(e[j + 1] - e[j]) for j in 1:n] ./ 2
    @assert abs(sum(w) - 1) < 1e-12 "lat weights sum to $(sum(w))"
    w
end

"""
    ocean_mean(r, yr) -> Float64 [m]

`<C·rsl> / <C>` at a stored epoch: relative sea level averaged over the ocean of
that moment. NOT the `bsl` diagnostic — `bsl` is the model's own eustatic
bookkeeping, this is an area average of the solved field.
"""
function ocean_mean(r::PubRun, yr::Real)
    C   = slice_at(r, "C_ocean", yr)
    rsl = slice_at(r, "rsl", yr)
    w   = lat_weights(r.lat)
    nlon = length(r.lon)
    sph(F) = sum((sum(@view F[:, j]) / nlon) * w[j] for j in eachindex(w))
    sph(C .* rsl) / sph(C)
end

"""
    polar_motion(r) -> (m1, m2) | nothing

The rotation vector perturbation on the full time axis. `nothing` for the VILMA1
runs, whose backend does not write it — a caller degrades rather than drawing an
empty line that would read as disagreement between the two codes.
"""
function polar_motion(r::PubRun)
    (haskey(r.ds, "rot_m_re") && haskey(r.ds, "rot_m_im")) || return nothing
    (Float64.(r.ds["rot_m_re"][:]), Float64.(r.ds["rot_m_im"][:]))
end

"""
    site_series(r, name, lat, lon) -> (values, dlat, dlon)

One site's full time series, with the offset of the grid cell from the requested
point. Sites are those of `Published.SITES` and are extracted at reduction; a
coordinate that is not on that list raises, because the underlying grid is not
in the published file.
"""
function site_series(r::PubRun, name::AbstractString, lat::Real, lon::Real)
    name == "rsl" || error("site_series: only rsl is published (asked for \"$name\")")
    slat = Float64.(r.ds["site_lat"][:]); slon = Float64.(r.ds["site_lon"][:])
    k = findfirst(i -> isapprox(slat[i], lat; atol = 1e-6) &&
                       isapprox(slon[i], lon; atol = 1e-6), eachindex(slat))
    k === nothing && error("""
        site ($(lat), $(lon)) is not in the published set.
        Published sites: $(join(String.(r.ds["site_name"][:]), ", ")).
        Add it to Published.SITES and re-run analysis/make_reduced.jl.
        """)
    (Float64.(r.ds["rsl_site"][k, :]),
     Float64(r.ds["site_dlat"][k]), Float64(r.ds["site_dlon"][k]))
end

"""
    regrid(field, lon_src, lat_src, lon_dst, lat_dst) -> Matrix

Bilinear interpolation of an oriented (lon, lat) field onto another oriented
grid. Longitude wraps; latitude clamps at the poles.

Kept here because the published files preserve each run's own Gauss grid, so a
reader comparing two truncations offline needs it for the same reason the
reduction did.
"""
function regrid(field::AbstractMatrix, lon_src, lat_src, lon_dst, lat_dst)
    out = similar(field, length(lon_dst), length(lat_dst))
    ns = length(lon_src)
    for (jd, la) in enumerate(lat_dst)
        j1 = searchsortedlast(lat_src, la)
        j1 = clamp(j1, 1, length(lat_src) - 1)
        tj = (la - lat_src[j1]) / (lat_src[j1 + 1] - lat_src[j1])
        tj = clamp(tj, 0.0, 1.0)
        for (id, lo) in enumerate(lon_dst)
            i1 = searchsortedlast(lon_src, lo)
            if i1 < 1 || i1 >= ns                # wrap across the −180/180 seam
                i1 = ns
                i2 = 1
                span = (lon_src[1] + 360) - lon_src[ns]
                ti = (lo - lon_src[ns] + (lo < lon_src[1] ? 360 : 0)) / span
            else
                i2 = i1 + 1
                ti = (lo - lon_src[i1]) / (lon_src[i2] - lon_src[i1])
            end
            ti = clamp(ti, 0.0, 1.0)
            out[id, jd] = (1 - ti) * (1 - tj) * field[i1, j1] +
                          ti       * (1 - tj) * field[i2, j1] +
                          (1 - ti) * tj       * field[i1, j1 + 1] +
                          ti       * tj       * field[i2, j1 + 1]
        end
    end
    out
end

regrid(field::AbstractMatrix, src::PubRun, dst::PubRun) =
    size(field) == (length(dst.lon), length(dst.lat)) ? field :
    regrid(field, src.lon, src.lat, dst.lon, dst.lat)

"""
    deg01(d, lon, lat) -> (c0, c1norm, residual, frac_var)

Area-weighted least-squares fit of degrees 0 and 1 — the constant and the three
direction cosines — to an oriented (lon, lat) field, and what is left after
removing them.

This is here rather than in a figure script because it is the measurement that
decides what the VILMA2−VILMA1 residual *is*. The residual was a degree-1
FRAME difference before it was anything else: VILMA2's historical `cf`
frame carries no degree 1 in rsl at all, while VILMA1 computes geocenter motion,
and on the disc benchmark that was 98-99.7 % of the disagreement. Running both
in the `cm` frame removed it. An offline reader re-deriving the residual needs
the same decomposition, so it travels with the data.
"""
function deg01(d::AbstractMatrix, lon, lat)
    LO = repeat(reshape(deg2rad.(lon), :, 1), 1, length(lat))
    LA = repeat(reshape(deg2rad.(lat), 1, :), length(lon), 1)
    x = cos.(LA) .* cos.(LO); y = cos.(LA) .* sin.(LO); z = sin.(LA); w = cos.(LA)
    A  = hcat(vec(ones(size(d))), vec(x), vec(y), vec(z))
    sw = sqrt.(vec(w))
    c  = (A .* sw) \ (vec(d) .* sw)
    fit = reshape(A * c, size(d))
    r   = d .- fit
    (c[1], sqrt(sum(c[2:4] .^ 2)), r, 1 - sum(w .* r .^ 2) / sum(w .* d .^ 2))
end

"Model year → ka for a label: −21000 → 21, and 0 → 0 (never −0)."
ka(t::Real) = -t / 1000 + 0.0

"""
    sh_parts(d, lon, lat; lmax = 8) -> (rms, parts)

Area-weighted rms of an oriented (lon, lat) field and its real spherical-harmonic
content to degree `lmax`, by projection with the Gauss quadrature weights:
`parts[(l, m)] = (rms, fit)`, where `rms` is the area-weighted rms of the
degree-l, order-|m| part (cos and sin together) and `fit` that part as a field.
"""
function sh_parts(d::AbstractMatrix, lon, lat; lmax::Int = 8)
    w   = lat_weights(lat)
    W   = repeat(reshape(w, 1, :), length(lon), 1) ./ length(lon)
    x   = sind.(Float64.(lat))
    φ   = deg2rad.(Float64.(lon))
    wrms(f) = sqrt(sum(W .* f .^ 2))
    # unnormalised associated Legendre P_l^m(x), l, m ≤ lmax
    P = Dict{Tuple{Int,Int},Vector{Float64}}()
    for m in 0:lmax
        pmm = ones(length(x)); fct = 1.0
        for _ in 1:m
            pmm .*= -fct .* sqrt.(max.(0.0, 1 .- x .^ 2)); fct += 2
        end
        P[(m, m)] = pmm
        m < lmax && (P[(m + 1, m)] = x .* (2m + 1) .* pmm)
        for l in m+2:lmax
            P[(l, m)] = (x .* (2l - 1) .* P[(l - 1, m)] .- (l + m - 1) .* P[(l - 2, m)]) ./ (l - m)
        end
    end
    parts = Dict{Tuple{Int,Int},Tuple{Float64,Matrix{Float64}}}()
    for l in 0:lmax, m in 0:l
        fit = zeros(size(d))
        for trig in (m == 0 ? (cos,) : (cos, sin))
            B = trig.(m .* φ) * P[(l, m)]'          # (lon, lat)
            nb = sum(W .* B .^ 2)
            nb > 0 && (fit .+= (sum(W .* d .* B) / nb) .* B)
        end
        parts[(l, m)] = (wrms(fit), fit)
    end
    (wrms(d), parts)
end

# ---------------------------------------------------------------------------
# The sea-level budget, precomputed
# ---------------------------------------------------------------------------
"""
    budget_series(r) -> NamedTuple

The deglacial sea-level budget over `BUDGET_ROWS` levels, all in barystatic
metres. Computed at reduction from the run's own reference ice and bed — see
`analysis/Budget.jl` for what each column means and why the decomposition is not
circular — and stored, because it needs fields at every level rather than at the
seven published epochs.
"""
budget_series(r::PubRun) = (
    time       = Float64.(r.ds["budget_time"][:]),
    rsl_mean   = Float64.(r.ds["budget_rsl_mean"][:]),
    coded      = Float64.(r.ds["budget_coded"][:]),
    required   = Float64.(r.ds["budget_required"][:]),
    target     = Float64.(r.ds["budget_target"][:]),
    flooding   = Float64.(r.ds["budget_flooding"][:]),
    ocean_frac = Float64.(r.ds["budget_ocean_frac"][:]),
    ref_ocean_frac = Float64(r.ds.attrib["budget_ref_ocean_frac"]),
)

"""
    budget_at(r) -> NamedTuple

The present-day row alone — the one the paper quotes.
"""
function budget_at(r::PubRun)
    b = budget_series(r)
    k = length(b.time)
    (; time = b.time[k], delivered = b.rsl_mean[k], target = b.target[k],
       closure = b.rsl_mean[k] / b.target[k],
       mask_term = b.coded[k] - b.required[k],
       required_plus_flooding = b.required[k] + b.flooding[k],
       coded_plus_flooding = b.coded[k] + b.flooding[k],
       ocean_frac = b.ocean_frac[k])
end

"""
    raw_series(r) -> NamedTuple

The three curves of the raw-field diagnostic, on the full time axis: peak |rsl|,
relative sea level averaged over the ocean of the moment, and the rms of the
field about that mean — the fingerprint with the mass term removed.

These use cos(lat) area weights rather than the Gaussian quadrature weights of
`ocean_mean`, which is the convention `v2v1_series_maps.jl` has always used; the
two answer the same question to within the difference between the weightings and
are deliberately not mixed.
"""
raw_series(r::PubRun) = (
    time    = r.time,
    max_abs = Float64.(r.ds["max_abs_rsl"][:]),
    omean   = Float64.(r.ds["ocean_mean_rsl"][:]),
    pattern = Float64.(r.ds["ocean_pattern_rms"][:]),
)

# ---------------------------------------------------------------------------
# The VILMA2 − VILMA1 residual, precomputed
# ---------------------------------------------------------------------------
"""
    load_v2v1() -> Vector{NamedTuple}

One entry per V2−V1 pair, each carrying the residual decomposition on
`V2V1_SERIES_TIMES` and the raw-field residual on the full time axis. Both need
two runs' fields differenced at every level, so neither can be rebuilt from the
published epochs.
"""
function load_v2v1()
    require_published(v2v1_path())
    ds = NCDataset(v2v1_path())
    try
        out = NamedTuple[]
        for k in eachindex(ds["pair_tag"][:])
            push!(out, (
                tag   = String(ds["pair_tag"][k]),
                label = String(ds["pair_label"][k]),
                v2_id = String(ds["pair_v2"][k]),
                v1_id = String(ds["pair_v1"][k]),
                time  = Float64.(ds["tv2v1_time"][:]),
                rms   = Float64.(ds["v2v1_rms"][k, :]),
                rel   = Float64.(ds["v2v1_rel"][k, :]),
                deg0  = Float64.(ds["v2v1_deg0"][k, :]),
                deg1  = Float64.(ds["v2v1_deg1"][k, :]),
                cmis  = Float64.(ds["v2v1_cmis"][k, :]),
                raw_time = Float64.(ds["time"][:]),
                raw_rms  = Float64.(ds["v2v1_rms_raw"][k, :]),
            ))
        end
        out
    finally
        close(ds)
    end
end

# ---------------------------------------------------------------------------
# Block E — the sensitivity axes, as numbers
# ---------------------------------------------------------------------------
"""
    load_blockE() -> (e1, e3, meta)

Block E is reduced to scalars: every quantity its two figures plot is a single
number per run, and the differences behind them (`rms` against the most-resolved
member, across four different Gauss grids) are computed at reduction. `e1` and
`e3` are vectors of NamedTuples in plotting order.
"""
function load_blockE()
    require_published(blockE_path())
    ds = NCDataset(blockE_path())
    try
        e1 = [(; id = String(ds["E1_id"][i]), lmax = Int(ds["E1_lmax"][i]),
                 ocean_mean = Float64(ds["E1_ocean_mean"][i]),
                 target     = Float64(ds["E1_target"][i]),
                 closure    = Float64(ds["E1_closure"][i]),
                 rms_rsl    = Float64(ds["E1_rms_rsl"][i]),
                 rms_upl    = Float64(ds["E1_rms_upl"][i]),
                 wall_min   = Float64(ds["E1_wall_min"][i]))
              for i in eachindex(ds["E1_lmax"][:])]
        e3 = [(; id = String(ds["E3_id"][i]), tol = Float64(ds["E3_tol"][i]),
                 label    = String(ds["E3_label"][i]),
                 ne3d       = Int(ds["E3_ne3d"][i]),
                 ocean_mean = Float64(ds["E3_ocean_mean"][i]),
                 rms_rsl  = Float64(ds["E3_rms_rsl"][i]),
                 rms_upl  = Float64(ds["E3_rms_upl"][i]),
                 wall_min = Float64(ds["E3_wall_min"][i]))
              for i in eachindex(ds["E3_tol"][:])]
        meta = (; ref_lmax = Int(ds.attrib["E1_reference_lmax"]),
                  ref_tol_id = String(ds.attrib["E3_reference_id"]),
                  ne_tot = Int(ds.attrib["E3_ne_total"]))
        (e1, e3, meta)
    finally
        close(ds)
    end
end

# ---------------------------------------------------------------------------
# Block A — the community benchmarks
# ---------------------------------------------------------------------------
"""
    BenchRun

One benchmark run, reduced. The driver's stdout is not published: what is
published is the numbers that were parsed out of it (`tables`) together with the
column tables the driver dumped alongside it (`dumps`), plus the verdict and the
SLURM accounting as metadata. A reader gets the measurements without having to
parse a DKRZ job epilogue.
"""
struct BenchRun
    id::String
    verdict::Union{String,Nothing}
    slurm::Dict{String,Any}
    dumps::Dict{String,Matrix{Float64}}
    tables::Dict{String,Any}
end

Base.show(io::IO, r::BenchRun) = print(io, "BenchRun($(r.id), $(something(r.verdict, "no verdict")))")

passed(r::BenchRun) = r.verdict !== nothing && startswith(r.verdict, "PASS")
verdict(r::BenchRun) = something(r.verdict, "(no verdict — job pending or crashed)")

"One of the driver's dumped column tables; a 0×0 matrix when that run has none."
dumptable(r::BenchRun, name::AbstractString) =
    get(r.dumps, replace(String(name), ".txt" => ""), Matrix{Float64}(undef, 0, 0))
hasdump(r::BenchRun, name::AbstractString) =
    haskey(r.dumps, replace(String(name), ".txt" => ""))

"""
    load_blockA() -> Vector{BenchRun}

Every reduced benchmark run, sorted by id.
"""
function load_blockA()
    require_published(blockA_path())
    ds = NCDataset(blockA_path())
    try
        out = BenchRun[]
        for id in sort(collect(keys(ds.group)))
            g = ds.group[id]
            v = haskey(g.attrib, "verdict") ? String(g.attrib["verdict"]) : nothing
            slurm = Dict{String,Any}()
            for k in ("jobid", "elapsed_s", "wall_s", "maxrss_kb")
                haskey(g.attrib, k) && (slurm[k] = Int(g.attrib[k]))
            end
            haskey(g.attrib, "node") && (slurm["node"] = String(g.attrib["node"]))
            dumps  = Dict{String,Matrix{Float64}}()
            tables = Dict{String,Any}()
            for vn in keys(g)
                startswith(vn, "dump_")  && (dumps[vn[6:end]]  = Float64.(g[vn][:, :]); continue)
                startswith(vn, "table_") && (tables[vn[7:end]] = _blockA_value(g[vn]))
            end
            push!(out, BenchRun(id, v, slurm, dumps, tables))
        end
        out
    finally
        close(ds)
    end
end

"A parsed table comes back as a matrix, a vector or a string vector, as written."
function _blockA_value(v)
    if eltype(v) <: AbstractString
        return String.(v[:])
    end
    ndims(v) == 2 ? Float64.(v[:, :]) : Float64.(v[:])
end

_tab(r::BenchRun, k, default) = get(r.tables, k, default)
const _EMPTY = Matrix{Float64}(undef, 0, 0)

"A1 — fluid-limit and elastic loading Love numbers, as the driver printed them."
love(r::BenchRun) = (; fluid = _tab(r, "love_fluid", _EMPTY),
                       elastic = _tab(r, "love_elastic", _EMPTY))

"A2 — the VE centre uplift table and the profile error summaries."
disc(r::BenchRun) = (; centre = _tab(r, "disc_centre", _EMPTY),
                       u_summary = _tab(r, "disc_u_summary", Float64[]),
                       n_summary = _tab(r, "disc_n_summary", Float64[]))

"A3 caseA — peak and max relative error per field."
martinec_caseA(r::BenchRun) = (; fields = _tab(r, "caseA_fields", String[]),
                                 peaks = _tab(r, "caseA_peaks", Float64[]),
                                 errs  = _tab(r, "caseA_errs", Float64[]))

"A3 C2/D3/E2/F1 — the peak-normalized max error table and the mass diagnostics."
sle(r::BenchRun) = (; figs = _tab(r, "sle_figs", String[]),
                      table = _tab(r, "sle_table", _EMPTY),
                      eustatic = _tab(r, "sle_eustatic", Float64[]),
                      ocean = _tab(r, "sle_ocean", Float64[]),
                      iters = _tab(r, "sle_iters", Float64[]))

"A4 — polar motion against Spada Table 14, per load geometry, and the tidal Love numbers."
function rotation(r::BenchRun)
    names = _tab(r, "rotation_bodies", String[])
    bodies = Dict{String,Matrix{Float64}}(b => _tab(r, "rotation_table_" * b, _EMPTY) for b in names)
    gnames = _tab(r, "rotation_G_bodies", String[])
    gvals  = _tab(r, "rotation_G", Float64[])
    G = Dict{String,Float64}(zip(gnames, gvals))
    kT = _tab(r, "rotation_kTe", [NaN]); ks = _tab(r, "rotation_ks", [NaN])
    (; bodies, G, kTe = isempty(kT) ? NaN : kT[1], k_s = isempty(ks) ? NaN : ks[1])
end

"A5 — the LVZ centre series and the 200 yr values against Weerdesteijn."
function lvz(r::BenchRun)
    names = _tab(r, "lvz_case_names", String[])
    vals  = _tab(r, "lvz_case_values", _EMPTY)
    cases = Dict{String,NTuple{3,Float64}}(
        n => (vals[i, 1], vals[i, 2], vals[i, 3]) for (i, n) in enumerate(names))
    (; t = _tab(r, "lvz_t", Float64[]), u = _tab(r, "lvz_u", Float64[]), cases,
       amplification = _tab(r, "lvz_amplification", Float64[]))
end

"""
    linear_nodes(y; rtol = 1e-9) -> Vector{Int}

Indices at which a piecewise-linear column carries independent data: the two
endpoints plus every interior point where the slope changes.

The Spada disc reference is distributed on a 0.1° θ grid but is piecewise linear
on a ~1° one, so only its breakpoints are data and comparing against every point
charges the model for the reference's own interpolation. Kept in the published
reader because the published dump tables carry that reference verbatim, so an
offline reader needs the same discrimination the paper's numbers rest on.
"""
function linear_nodes(y::AbstractVector; rtol = 1e-9)
    n = length(y)
    n < 3 && return collect(1:n)
    curv = abs.(diff(diff(y)))
    m = maximum(curv)
    m == 0 && return [1, n]
    brk = findall(>(rtol * m), curv) .+ 1
    sort!(unique(vcat(1, brk, n)))
end

# ---------------------------------------------------------------------------
# Block B — the lateral-viscosity benchmark (Klemann et al. setup note)
# ---------------------------------------------------------------------------
# Published as the driver's own protocol files, unchanged: they are small text
# tables in the format the benchmark asks participants to submit, so a copy is
# both the reduced form and the submission. The shape below is what the reducer
# checks before it copies anything.

"Tests A-i, B-i, C-i and the two load histories."
const BENCH_TESTS    = ["A", "B", "C"]
const BENCH_FORCINGS = ["heav", "ramp"]
"The protocol's output epochs [kyr] (tests/bench_visc3d.f90 `T_OUT_KYR`)."
const BENCH_EPOCHS_KYR = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 11.0, 12.0, 15.0, 20.0, 50.0, 100.0]
"Points on the cross-section, structure centre to antipode every 0.1°."
const BENCH_NPT = 1801
"The submission's truncation; its files sit at the top of `bench_visc3d/`."
const BENCH_LMAX = 128
"The convergence truncations, run for the axisymmetric tests only, under `convergence/l<L>/`."
const BENCH_CONV_LMAX  = [64, 256]
const BENCH_CONV_TESTS = ["A", "B"]

"Directory holding the protocol files of truncation `lmax`."
bench_case_dir(lmax) = lmax == BENCH_LMAX ? blockB_dir() :
                       joinpath(blockB_dir(), "convergence", "l$lmax")

"""
The normal-mode Love numbers of M3-L70-V01 (giapy; Spada et al. 2011), copied from
the model repo at reduction: the exact reference for test A-i (analysis/NormalMode.jl).
"""
bench_love_path() = joinpath(blockB_dir(), "mod_M3-L70-V01")

bench_disp_name(t, f) = "disp_VILMA2_$(t)-i_$(f).txt"
bench_gcm_name(t, f)  = "gcm_VILMA2_$(t)-i_$(f).txt"
bench_cost_name(t, f) = "cost_$(t)-i_$(f).txt"

"""
    bench_cost(t, f; lmax) -> (wall_s, nseg, threads) | nothing

Wall time of the whole chain of segment jobs of one case, from the jobs' own
accounting, written at reduction.
"""
function bench_cost(t, f; lmax = BENCH_LMAX)
    p = joinpath(bench_case_dir(lmax), bench_cost_name(t, f))
    isfile(p) || return nothing
    kv = Dict(split(l)[1] => parse(Float64, split(l)[2]) for l in eachline(p) if length(split(l)) == 2)
    (wall_s = kv["wall_s"], nseg = Int(kv["segments"]), threads = Int(get(kv, "threads", 0)))
end

"""
One benchmark case. `disp` columns: lon, lat, distance along the section [deg],
time [kyr], u_r, u_ϑ, u_φ [m], δφ [m² s⁻²]. `gcm` columns: time [kyr], u_x, u_y,
u_z [m], the geocentre motion u_CF − u_CM. `lmax` is read from the header.
"""
struct BenchV3D
    test::String
    forcing::String
    lmax::Int
    disp::Matrix{Float64}
    gcm::Matrix{Float64}
end

"Numeric rows of a protocol file (whitespace-delimited, `#` header lines)."
function read_bench_table(path)
    rows = [parse.(Float64, split(strip(l))) for l in eachline(path)
            if !startswith(l, "#") && !isempty(strip(l))]
    reduce(vcat, permutedims.(rows))
end

"`lmax` from the protocol header line `# C-i (heav), lmax=128`."
function bench_lmax(path)
    for l in eachline(path)
        m = match(r"lmax=(\d+)", l)
        m === nothing || return parse(Int, m[1])
        startswith(l, "#") || break
    end
    error("$path: no lmax in the header")
end

"""
    load_blockB(lmax = BENCH_LMAX) -> Dict{Tuple{String,String},BenchV3D}

Every published case at truncation `lmax`, keyed `(test, forcing)`. A case that
was not published is absent rather than empty.
"""
function load_blockB(lmax = BENCH_LMAX)
    dir = bench_case_dir(lmax)
    require_published(dir)
    out = Dict{Tuple{String,String},BenchV3D}()
    for t in BENCH_TESTS, f in BENCH_FORCINGS
        pd = joinpath(dir, bench_disp_name(t, f))
        pg = joinpath(dir, bench_gcm_name(t, f))
        (isfile(pd) && isfile(pg)) || continue
        out[(t, f)] = BenchV3D(t, f, bench_lmax(pd), read_bench_table(pd), read_bench_table(pg))
    end
    out
end

end # module
