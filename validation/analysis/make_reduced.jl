#!/usr/bin/env julia
#
# make_reduced.jl — distil `runs/` into `data/published/`, the small set of
# files the figures and tables are made from.
#
#   julia --project=. analysis/make_reduced.jl [BLOCK ...]
#
# with BLOCK one or more of D, V2V1, E, A, B; no argument does all of them. A
# block whose runs are missing or incomplete is skipped.
#
# This is the only script that reads `runs/`. Everything downstream goes through
# `analysis/Published.jl` against `data/published/`.
#
#   * Subset: 2-D fields are kept at `Published.EPOCHS` only.
#   * Derive: `z_bed` is not stored (it is `z_bed_eq - rsl`, asserted below).
#   * Precompute: whatever sweeps every time level is stored as a series: the
#     sea-level budget, the V2−V1 residual decomposition, the raw-field
#     diagnostics, the per-site histories, and all of block E.
#   * Parse: block A exists as driver stdout and dump tables; the numbers parsed
#     out of them are stored.
#
# Fields are written the way they are plotted: longitude −180…180, latitude
# south to north. A full pass reads every time level of the block D files
# twice; allow ten to twenty minutes and run it in a job (./run.sh submit reduce).

using NCDatasets
using Printf
using Statistics

include(joinpath(@__DIR__, "BlockD.jl"))
include(joinpath(@__DIR__, "BlockA.jl"))
include(joinpath(@__DIR__, "Budget.jl"))
include(joinpath(@__DIR__, "Published.jl"))
using .BlockD
using .BlockA
using .Budget
# NOT `using .Published`: it exports the same names as BlockD, Budget and BlockA
# by design -- it is the same interface served from the reduced files instead of
# from runs/. Bringing both into scope would make every one of those names
# ambiguous. Here the run-side modules are the ones in use, and Published is
# referenced qualified, for the schema constants and `deg01`.

const REPO = dirname(@__DIR__)
const RUNS = joinpath(REPO, "runs")
const OUT  = Published.published_dir()

# NetCDF4 with shuffle+deflate. The published fields are smooth (rsl), mostly
# zero (h_ice) or binary (C_ocean), so this is the difference between a set that
# can be attached to a paper and one that cannot.
const DEFLATE = (deflatelevel = 4, shuffle = true)

"""
Fields are published at the precision the model wrote them, not at Float32.

Float32 looks like an obvious halving and is not: `z_bed` is derived downstream
as `z_bed_eq - rsl`, and at bedrock elevations of a few thousand metres Float32
rounding is 2^-11 m. That is invisible in a map and NOT invisible in the tables,
where the V2-V1 uplift residual is quoted to 1e-4 mm/yr and is itself a difference
of two rsl fields 500 yr apart. Storing the primary fields exactly means the
published data reproduces the paper's numbers digit for digit rather than nearly,
and deflate recovers most of the difference anyway because these fields are
smooth, mostly zero, or binary.
"""
const FIELD_T = Float64

say(args...) = (println("\n== ", args...); flush(stdout))

# The block D production matrix: two Earth structures, two solvers.
const BLOCKD_RUNS = [
    (id = "D-V2-3D.l170.fe",  root = "blockD_v2",    solver = "VILMA2"),
    (id = "D-V2-1Db.l170.fe", root = "blockD_v2",    solver = "VILMA2"),
    (id = "D-V1-3D.l170.fe",  root = "blockD_v1", solver = "VILMA1"),
    (id = "D-V1-1Db.l170.fe", root = "blockD_v1", solver = "VILMA1"),
    # V7: the 3-D run with the toroidal field on; D-V2-3D.l170.fe is its control.
    (id = "D-V2-3D.l170.tor",   root = "blockD_tor", solver = "VILMA2"),
]

# The pruned V2−V1 list: every pair on VILMA2's own LGM equilibrium
# (i_eq = 2) with the degree-1 frames matched (deg1_frame = "cm"). See
# analysis/check_v2v1.jl for why those two choices and no others. 3d_tor is V7:
# the 3-D pair again, with VILMA2's toroidal coupling on.
const V2V1_PAIRS = [
    (tag = "3d",  v2 = "D-V2-3D.l170.fe",  v1 = "D-V1-3D.l170.fe",
     label = "3-D Bagge (2021) viscosity"),
    (tag = "1db", v2 = "D-V2-1Db.l170.fe", v1 = "D-V1-1Db.l170.fe",
     label = "1Db, the radial mean of that field"),
    (tag = "3d_tor", v2 = "D-V2-3D.l170.tor", v1 = "D-V1-3D.l170.fe",
     label = "3-D, VILMA2 with the toroidal coupling"),
]

# Block E, as run_blockE.jl declares it: each axis measured against its OWN
# most-resolved member, not against the production setting.
const E1 = [
    (lmax =  64, root = "blockE",    id = "D-V2-3D.l64.cold"),
    (lmax = 128, root = "blockE",    id = "D-V2-3D.l128.cold"),
    (lmax = 170, root = "blockD_v2", id = "D-V2-3D.l170.fe"),
    (lmax = 256, root = "blockE",    id = "D-V2-3D.l256.cold"),
]
const E3 = [
    (tol = 0.0,    label = "0 (all)", root = "blockE", id = "E3-tolall.l170"),
    (tol = 1.0e-3, label = "10⁻³",    root = "blockE", id = "E3-tol1e-3.l170"),
    (tol = 1.0e-2, label = "10⁻²",    root = "blockE", id = "E3-tol1e-2.l170"),
    (tol = 1.0e-1, label = "10⁻¹",    root = "blockE", id = "E3-tol1e-1.l170"),
]

complete(r) = r !== nothing && length(r.time) >= 3 && r.time[end] >= -1.0
rundir(root, id) = joinpath(RUNS, root, id)

"Wall clock [s] of a run directory, from the accounting lines at the end of out.out."
function wall_s(root, id)
    f = joinpath(RUNS, root, id, "out.out")
    isfile(f) ? get(BlockA.slurm_meta(readlines(f)), "elapsed_s", nothing) : nothing
end

"Write a variable with a fresh, uniquely named dimension per axis."
function putvar(g, name, data; attrib = Dict{String,Any}(), kwargs...)
    dims = String[]
    for (k, n) in enumerate(size(data))
        d = "$(name)_d$(k)"
        haskey(g.dim, d) || defDim(g, d, n)
        push!(dims, d)
    end
    v = defVar(g, name, eltype(data), Tuple(dims); attrib, kwargs...)
    v[:] = data
    v
end

# =============================================================================
# Block D — the deglaciation runs
# =============================================================================
"""
Area weights on the oriented grid, in the cos(lat) convention that the raw-field
diagnostics have always used. Deliberately NOT the Gaussian quadrature weights
of `ocean_mean`: the two answer the same question to within the difference
between the weightings, and mixing them would silently change a published curve.
"""
cos_wts(r) = repeat(reshape(cosd.(r.lat), 1, :), length(r.lon), 1)

function reduce_blockD()
    say("block D — the deglaciation runs")
    mkpath(joinpath(OUT, "blockD"))
    done = String[]
    for c in BLOCKD_RUNS
        dir = rundir(c.root, c.id)
        r = load_deglac(dir; id = c.id)
        if !complete(r)
            @warn "skipping $(c.id) (missing or incomplete under $dir)"
            continue
        end
        @info "reducing $(c.id): $(length(r.time)) levels, $(length(r.lon))×$(length(r.lat))"
        reduce_one_blockD(r, c)
        push!(done, c.id)
    end
    isempty(done) && @warn "no block D runs were reduced"
    done
end

function reduce_one_blockD(r, c)
    path = joinpath(OUT, "blockD", c.id * ".nc")
    nlon, nlat, nt = length(r.lon), length(r.lat), length(r.time)

    # --- the epoch subset -------------------------------------------------
    # Snap each requested epoch onto the run's own time axis and refuse to
    # publish a field that is not the epoch it claims to be.
    eidx = Int[]
    for yr in Published.EPOCHS
        i = nearest_t(r, yr)
        abs(r.time[i] - yr) <= Published.EPOCH_TOL ||
            error("$(c.id): nearest level to $yr yr is $(r.time[i]) yr, further than " *
                  "$(Published.EPOCH_TOL) yr. This run's axis cannot carry Published.EPOCHS.")
        push!(eidx, i)
    end
    epoch = r.time[eidx]

    rsl   = Array{FIELD_T}(undef, nlon, nlat, length(eidx))
    hice  = similar(rsl)
    cocn  = similar(rsl)
    for (k, i) in enumerate(eidx)
        rsl[:, :, k]  = slice_at(r, "rsl", r.time[i])
        hice[:, :, k] = slice_at(r, "h_ice", r.time[i])
        cocn[:, :, k] = slice_at(r, "C_ocean", r.time[i])
    end
    zbeq = FIELD_T.(orient_eq(r, "z_bed_eq"))
    hieq = FIELD_T.(orient_eq(r, "h_ice_eq"))

    # z_bed is dropped and derived downstream. Assert the identity on every
    # published epoch rather than trust the variable's long_name -- it is a
    # quarter of the published volume riding on one documented relation, and the
    # check has to be exact, which is the other half of why FIELD_T is Float64.
    for (k, i) in enumerate(eidx)
        zb = slice_at(r, "z_bed", r.time[i])
        d = maximum(abs.(zbeq .- rsl[:, :, k] .- zb))
        d == 0 || error("$(c.id): z_bed != z_bed_eq - rsl at $(r.time[i]) yr " *
                        "(max |diff| = $d). Published.slice_at derives z_bed and " *
                        "would be wrong.")
    end

    # --- the series on the full axis --------------------------------------
    # One pass over every level: peak |rsl|, ocean-mean rsl and the rms about
    # that mean. These cannot be rebuilt from seven epochs, and they are what
    # the raw-field diagnostic plots.
    mx = zeros(nt); om = zeros(nt); pr = zeros(nt)
    w0 = cos_wts(r)
    for i in 1:nt
        f = slice_at(r, "rsl", r.time[i])
        cc = slice_at(r, "C_ocean", r.time[i])
        w = w0 .* (cc .> 0.5)
        sw = max(sum(w), eps())
        m = sum(w .* f) / sw
        mx[i] = maximum(abs.(f))
        om[i] = m
        pr[i] = sqrt(sum(w .* (f .- m) .^ 2) / sw)
        (i % 50 == 0 || i == nt) && (@printf("    series %3d/%d\r", i, nt); flush(stdout))
    end
    println()

    sites = [site_series(r, "rsl", lat, lon) for (_, lat, lon) in Published.SITES]
    rsl_site = permutedims(reduce(hcat, [s[1] for s in sites]))   # (site, time)

    bud = budget_series(joinpath(r.dir, "out.nc"); nrow = Published.BUDGET_ROWS)
    pm  = polar_motion(r)

    # --- write ------------------------------------------------------------
    isfile(path) && rm(path)
    ds = NCDataset(path, "c")
    try
        defDim(ds, "lon", nlon); defDim(ds, "lat", nlat)
        defDim(ds, "epoch", length(eidx)); defDim(ds, "time", nt)
        defDim(ds, "site", length(Published.SITES))
        defDim(ds, "tbudget", length(bud.time))

        ds.attrib["title"]   = "VILMA2 paper — reduced deglaciation output, $(c.id)"
        ds.attrib["id"]      = c.id
        ds.attrib["solver"]  = c.solver
        ds.attrib["visc"]    = r.visc
        ds.attrib["lmax"]    = r.lmax
        ds.attrib["source_run"] = relpath(r.dir, REPO)
        ds.attrib["source_file_bytes"] = filesize(joinpath(r.dir, "out.nc"))
        ds.attrib["generated_by"] = "analysis/make_reduced.jl"
        ds.attrib["generated_on"] = Dates_now()
        ds.attrib["Conventions"]  = "CF-1.8"
        ds.attrib["comment"] = "Reduced from the full model output: 2-D fields at " *
            "selected epochs only; z_bed omitted and equal to z_bed_eq - rsl; " *
            "diagnostics that need every time level are precomputed as series."
        for k in ("ne", "ne3d", "visc3d_tol", "l_visc_3d")
            haskey(r.ds.attrib, k) && (ds.attrib[k] = r.ds.attrib[k])
        end
        w = wall_s(c.root, c.id); w === nothing || (ds.attrib["elapsed_s"] = w)
        ds.attrib["budget_ref_ocean_frac"] = bud.ref_ocean_frac

        v = defVar(ds, "lon", Float64, ("lon",);
                   attrib = ["units" => "degrees_east", "long_name" => "longitude",
                             "axis" => "X"]); v[:] = r.lon
        v = defVar(ds, "lat", Float64, ("lat",);
                   attrib = ["units" => "degrees_north", "long_name" => "latitude",
                             "axis" => "Y"]); v[:] = r.lat
        v = defVar(ds, "epoch", Float64, ("epoch",);
                   attrib = ["units" => "years", "long_name" =>
                             "epoch of the published 2-D fields (negative = before present)"])
        v[:] = epoch
        v = defVar(ds, "time", Float64, ("time",);
                   attrib = ["units" => "years", "axis" => "T", "long_name" =>
                             "coupling time axis of the published series"]); v[:] = r.time

        fld(name, data, units, long) = begin
            vv = defVar(ds, name, FIELD_T, ("lon", "lat", "epoch");
                        attrib = ["units" => units, "long_name" => long], DEFLATE...)
            vv[:, :, :] = data
        end
        fld("rsl", rsl, "m", "Relative sea-level change (full field, vs the run's reference state)")
        fld("h_ice", hice, "m", "Grounded-ice thickness")
        fld("C_ocean", cocn, "1", "Ocean function (1 ocean / 0 land)")

        v = defVar(ds, "z_bed_eq", FIELD_T, ("lon", "lat");
                   attrib = ["units" => "m",
                             "long_name" => "Reference (equilibrium) bedrock elevation"],
                   DEFLATE...); v[:, :] = zbeq
        v = defVar(ds, "h_ice_eq", FIELD_T, ("lon", "lat");
                   attrib = ["units" => "m",
                             "long_name" => "Reference (equilibrium) grounded-ice thickness"],
                   DEFLATE...); v[:, :] = hieq

        ser(name, data, units, long) = begin
            vv = defVar(ds, name, Float64, ("time",);
                        attrib = ["units" => units, "long_name" => long]); vv[:] = data
        end
        ser("bsl", Float64.(r.ds["bsl"][:]), "m",
            "Barystatic sea level vs reference (the model's own eustatic bookkeeping)")
        if pm !== nothing
            ser("rot_m_re", pm[1], "1", "Polar motion m1 (rotation, real part)")
            ser("rot_m_im", pm[2], "1", "Polar motion m2 (rotation, imaginary part)")
        end
        ser("max_abs_rsl", mx, "m", "Peak |rsl| over the globe")
        ser("ocean_mean_rsl", om, "m",
            "rsl averaged over the ocean of the moment (cos-lat weights)")
        ser("ocean_pattern_rms", pr, "m",
            "rms of rsl about its own ocean mean (the fingerprint, mass term removed)")

        v = defVar(ds, "site_name", String, ("site",);
                   attrib = ["long_name" => "site label"])
        v[:] = [String(n) for (n, _, _) in Published.SITES]
        for (nm, data, units, long) in (
                ("site_lat",  [lat for (_, lat, _) in Published.SITES], "degrees_north", "requested site latitude"),
                ("site_lon",  [lon for (_, _, lon) in Published.SITES], "degrees_east",  "requested site longitude"),
                ("site_dlat", [s[2] for s in sites], "degrees", "offset of the grid cell used, in latitude"),
                ("site_dlon", [s[3] for s in sites], "degrees", "offset of the grid cell used, in longitude"))
            vv = defVar(ds, nm, Float64, ("site",);
                        attrib = ["units" => units, "long_name" => long]); vv[:] = data
        end
        v = defVar(ds, "rsl_site", Float64, ("site", "time");
                   attrib = ["units" => "m", "long_name" =>
                             "Relative sea level at the nearest grid cell to each site"])
        v[:, :] = rsl_site

        bud_v(name, data, units, long) = begin
            vv = defVar(ds, name, Float64, ("tbudget",);
                        attrib = ["units" => units, "long_name" => long]); vv[:] = data
        end
        bud_v("budget_time", bud.time, "years", "time levels of the sea-level budget")
        bud_v("budget_rsl_mean", bud.rsl_mean, "m", "delivered ocean-mean relative sea level")
        bud_v("budget_coded", bud.coded, "m", "melt source as vilma_sle.f90 computes it")
        bud_v("budget_required", bud.required, "m",
              "the same source with each endpoint masked by its own coastline")
        bud_v("budget_target", bud.target, "m",
              "above-flotation volume of the reference ice lost, in barystatic metres")
        bud_v("budget_flooding", bud.flooding, "m",
              "subgrid term charging the ocean for water filling a newly flooded cell")
        bud_v("budget_ocean_frac", bud.ocean_frac, "1", "ocean fraction at that moment")
    finally
        close(ds)
    end
    @info @sprintf("  wrote %s (%.1f MB from %.0f MB)",
                   relpath(path, REPO), filesize(path) / 1e6,
                   filesize(joinpath(r.dir, "out.nc")) / 1e6)
    close(r.ds)
end

"A time-independent 2-D field, oriented the way the published files store them."
function orient_eq(r, name)
    a = Float64.(r.ds[name][:, :])
    reverse(circshift(a, (r.roll, 0)), dims = 2)
end

Dates_now() = Libc.strftime("%Y-%m-%d %H:%M:%S %Z", time())

# =============================================================================
# V2V1 — the VILMA2 − VILMA1 residual
# =============================================================================
function reduce_v2v1()
    say("V2−V1 — the VILMA2 against VILMA1 residual")
    pairs = NamedTuple[]
    for p in V2V1_PAIRS
        root2 = findfirst(c -> c.id == p.v2, BLOCKD_RUNS)
        root1 = findfirst(c -> c.id == p.v1, BLOCKD_RUNS)
        f = load_deglac(rundir(BLOCKD_RUNS[root2].root, p.v2); id = "V2.$(p.tag)")
        v = load_deglac(rundir(BLOCKD_RUNS[root1].root, p.v1); id = "V1.$(p.tag)")
        if !(complete(f) && complete(v))
            @warn "pair $(p.tag) incomplete — skipping"
            continue
        end
        @info "pair $(p.tag): V2 $(length(f.time)) steps, V1 $(length(v.time)) steps"

        # The decomposition, at the epochs the summary figure is drawn on.
        n = length(Published.V2V1_SERIES_TIMES)
        s = (time = zeros(n), rms = zeros(n), rel = zeros(n),
             deg0 = zeros(n), deg1 = zeros(n), cmis = zeros(n))
        for (k, t) in enumerate(Published.V2V1_SERIES_TIMES)
            d  = rsl_anomaly(f, t) .- rsl_anomaly(v, t)
            sg = maximum(abs.(rsl_anomaly(f, t)))
            c0, c1, _, _ = Published.deg01(d, f.lon, f.lat)
            cf = slice_at(f, "C_ocean", t); cv = slice_at(v, "C_ocean", t)
            s.time[k] = f.time[nearest_t(f, t)]
            s.rms[k]  = sqrt(mean(d .^ 2))
            s.rel[k]  = sg == 0 ? NaN : sqrt(mean(d .^ 2)) / sg
            s.deg0[k] = c0
            s.deg1[k] = c1
            s.cmis[k] = 100 * count(abs.(cf .- cv) .> 0.5) / length(cf)
            @printf("    decomposition %2d/%d\r", k, n); flush(stdout)
        end
        println()

        # The residual on the RAW field, every level. This is the diagnostic
        # that tells the opposite story from the anomaly — against the
        # present-day datum the residual appears to decay through the
        # deglaciation, on the raw field it grows monotonically, because the
        # datum itself carries most of the disagreement.
        nt = length(f.time)
        raw = zeros(nt)
        for i in 1:nt
            t = f.time[i]
            raw[i] = sqrt(mean((slice_at(f, "rsl", t) .- slice_at(v, "rsl", t)) .^ 2))
            (i % 50 == 0 || i == nt) && (@printf("    raw residual %3d/%d\r", i, nt); flush(stdout))
        end
        println()

        push!(pairs, (; p.tag, p.label, v2_id = p.v2, v1_id = p.v1, s, raw, time = f.time))
        close(f.ds); close(v.ds)
    end
    isempty(pairs) && (@warn "no V2−V1 pairs were reduced"; return)

    path = Published.v2v1_path()
    mkpath(dirname(path)); isfile(path) && rm(path)
    ds = NCDataset(path, "c")
    try
        np = length(pairs); ntf = length(Published.V2V1_SERIES_TIMES)
        nt = length(pairs[1].time)
        defDim(ds, "pair", np); defDim(ds, "tv2v1", ntf); defDim(ds, "time", nt)
        ds.attrib["title"] = "VILMA2 paper — VILMA2 minus VILMA1 residual"
        ds.attrib["generated_by"] = "analysis/make_reduced.jl"
        ds.attrib["generated_on"] = Dates_now()
        ds.attrib["comment"] = "Pairwise quantities: each entry needs both runs' full " *
            "fields differenced at every level, so none of it can be rebuilt from the " *
            "epochs published in blockD/. Sea level is referenced to each run's own " *
            "present day, so the residual at 0 ka is zero by construction."

        sv(name, vals; long = "") = begin
            v = defVar(ds, name, String, ("pair",); attrib = ["long_name" => long])
            v[:] = vals
        end
        sv("pair_tag",   [String(p.tag) for p in pairs]; long = "pair identifier")
        sv("pair_label", [String(p.label) for p in pairs]; long = "Earth structure of the pair")
        sv("pair_v2",    [String(p.v2_id) for p in pairs]; long = "VILMA2 run id")
        sv("pair_v1",    [String(p.v1_id) for p in pairs]; long = "VILMA1 run id")

        v = defVar(ds, "tv2v1_time", Float64, ("tv2v1",);
                   attrib = ["units" => "years",
                             "long_name" => "epochs of the residual decomposition"])
        v[:] = pairs[1].s.time
        v = defVar(ds, "time", Float64, ("time",);
                   attrib = ["units" => "years", "axis" => "T",
                             "long_name" => "coupling time axis of the raw residual"])
        v[:] = pairs[1].time

        pv(name, get, units, long) = begin
            vv = defVar(ds, name, Float64, ("pair", "tv2v1");
                        attrib = ["units" => units, "long_name" => long])
            for (k, p) in enumerate(pairs); vv[k, :] = get(p); end
        end
        pv("v2v1_rms",  p -> p.s.rms,  "m", "rms of (V2 − V1) in rsl vs present day")
        pv("v2v1_rel",  p -> p.s.rel,  "1", "that rms as a fraction of max |rsl| in V2")
        pv("v2v1_deg0", p -> p.s.deg0, "m", "degree-0 content of the residual")
        pv("v2v1_deg1", p -> p.s.deg1, "m", "degree-1 magnitude of the residual")
        pv("v2v1_cmis", p -> p.s.cmis, "%", "grid cells whose ocean functions disagree")

        vv = defVar(ds, "v2v1_rms_raw", Float64, ("pair", "time");
                    attrib = ["units" => "m", "long_name" =>
                              "rms of (V2 − V1) in the RAW rsl field, not the present-day anomaly"])
        for (k, p) in enumerate(pairs); vv[k, :] = p.raw; end
    finally
        close(ds)
    end
    @info @sprintf("  wrote %s (%.2f MB)", relpath(path, REPO), filesize(path) / 1e6)
end

# =============================================================================
# Block E — the sensitivity axes, reduced to numbers
# =============================================================================
"""
    against(run, ref) -> (rms_rsl_21ka, rms_uplift)

Difference of `run` from `ref`, both brought onto whichever of the two grids is
coarser, in rms over the global grid. Interpolating the coarse field UP instead
would invent structure it does not have and report the invention as convergence
error.
"""
function against(run, ref)
    coarse = length(run.lon) <= length(ref.lon) ? run : ref
    a  = BlockD.regrid(rsl_anomaly(run, -21000.0), run, coarse)
    b  = BlockD.regrid(rsl_anomaly(ref, -21000.0), ref, coarse)
    ua = BlockD.regrid(uplift_rate(run), run, coarse)
    ub = BlockD.regrid(uplift_rate(ref), ref, coarse)
    (sqrt(mean((a .- b) .^ 2)), sqrt(mean((ua .- ub) .^ 2)))
end

ne_attr(r, k) = haskey(r.ds.attrib, k) ? Int(r.ds.attrib[k]) : 0

function reduce_blockE()
    say("block E — the sensitivity axes")

    load(e) = begin
        r = load_deglac(rundir(e.root, e.id); id = e.id)
        complete(r) ? (; e..., run = r) : (@warn "block E: skipping $(e.id)"; nothing)
    end
    e1 = filter(!isnothing, [load(e) for e in E1])
    e3 = filter(!isnothing, [load(e) for e in E3])
    if length(e1) < 2 && length(e3) < 2
        @warn "block E has too few runs to reduce"
        return
    end

    rows1 = NamedTuple[]
    if length(e1) >= 2
        ref = e1[argmax([x.lmax for x in e1])]
        @info "E1: reference is ℓmax $(ref.lmax)"
        for x in e1
            d = x.lmax == ref.lmax ? (NaN, NaN) : against(x.run, ref.run)
            b = budget_at(joinpath(rundir(x.root, x.id), "out.nc"))
            push!(rows1, (; id = x.id, lmax = x.lmax,
                            ocean_mean = ocean_mean(x.run, 0.0),
                            target = b.target, closure = 100b.closure,
                            rms_rsl = d[1], rms_upl = d[2],
                            wall_min = something(wall_s(x.root, x.id), NaN) / 60))
            @info @sprintf("  ℓmax %3d: delivered %.2f m, target %.2f m, rms %.3f m",
                           x.lmax, rows1[end].ocean_mean, rows1[end].target, d[1])
        end
        refl = ref.lmax
    else
        refl = 0
    end

    rows3 = NamedTuple[]
    net = 0
    if length(e3) >= 2
        ref = e3[argmin([x.tol for x in e3])]
        @info "E3: reference is tol $(ref.tol) ($(ref.id))"
        net = ne_attr(e3[1].run, "ne")
        for x in e3
            d = x.id == ref.id ? (NaN, NaN) : against(x.run, ref.run)
            push!(rows3, (; id = x.id, tol = x.tol, label = x.label,
                            ne3d = ne_attr(x.run, "ne3d"),
                            ocean_mean = ocean_mean(x.run, 0.0),
                            rms_rsl = d[1], rms_upl = d[2],
                            wall_min = something(wall_s(x.root, x.id), NaN) / 60))
            @info @sprintf("  tol %-8s: ne3d %3d, rms %.4f m", x.label, rows3[end].ne3d, d[1])
        end
        ref3 = ref.id
    else
        ref3 = ""
    end

    path = Published.blockE_path()
    mkpath(dirname(path)); isfile(path) && rm(path)
    ds = NCDataset(path, "c")
    try
        ds.attrib["title"] = "VILMA2 paper — block E sensitivity, reduced to scalars"
        ds.attrib["generated_by"] = "analysis/make_reduced.jl"
        ds.attrib["generated_on"] = Dates_now()
        ds.attrib["comment"] = "Every quantity block E's figures plot is one number per " *
            "run. The rms differences are computed against each axis' own most-resolved " *
            "member (ℓmax 256 for E1, tol = 0 for E3) on the coarser grid of each pair, " *
            "which needs the full fields and therefore happens at reduction."
        ds.attrib["E1_reference_lmax"] = refl
        ds.attrib["E3_reference_id"]   = ref3
        ds.attrib["E3_ne_total"]       = net

        if !isempty(rows1)
            defDim(ds, "e1", length(rows1))
            v = defVar(ds, "E1_id", String, ("e1",)); v[:] = [String(r.id) for r in rows1]
            e1v(n, g, u, l) = (vv = defVar(ds, n, Float64, ("e1",);
                    attrib = ["units" => u, "long_name" => l]); vv[:] = [g(r) for r in rows1])
            e1v("E1_lmax", r -> Float64(r.lmax), "1", "spectral truncation")
            e1v("E1_ocean_mean", r -> r.ocean_mean, "m", "delivered ocean-mean rsl at present day")
            e1v("E1_target", r -> r.target, "m", "above-flotation volume of reference ice lost")
            e1v("E1_closure", r -> r.closure, "%", "delivered / target")
            e1v("E1_rms_rsl", r -> r.rms_rsl, "m", "rms difference from the reference truncation, rsl at 21 ka")
            e1v("E1_rms_upl", r -> r.rms_upl, "mm/yr", "rms difference from the reference truncation, uplift rate")
            e1v("E1_wall_min", r -> r.wall_min, "min", "wall clock for the full deglaciation")
        end
        if !isempty(rows3)
            defDim(ds, "e3", length(rows3))
            v = defVar(ds, "E3_id", String, ("e3",)); v[:] = [String(r.id) for r in rows3]
            v = defVar(ds, "E3_label", String, ("e3",)); v[:] = [String(r.label) for r in rows3]
            e3v(n, g, u, l) = (vv = defVar(ds, n, Float64, ("e3",);
                    attrib = ["units" => u, "long_name" => l]); vv[:] = [g(r) for r in rows3])
            e3v("E3_tol", r -> r.tol, "1", "visc3d_tol, in dex of lateral viscosity spread")
            e3v("E3_ne3d", r -> Float64(r.ne3d), "1", "radial elements advanced as laterally 3-D")
            e3v("E3_ocean_mean", r -> r.ocean_mean, "m", "delivered ocean-mean rsl at present day")
            e3v("E3_rms_rsl", r -> r.rms_rsl, "m", "rms difference from tol = 0, rsl at 21 ka")
            e3v("E3_rms_upl", r -> r.rms_upl, "mm/yr", "rms difference from tol = 0, uplift rate")
            e3v("E3_wall_min", r -> r.wall_min, "min", "wall clock for the full deglaciation")
        end
    finally
        close(ds)
    end
    for x in vcat(e1, e3); close(x.run.ds); end
    @info @sprintf("  wrote %s (%.2f MB)", relpath(path, REPO), filesize(path) / 1e6)
end

# =============================================================================
# Block A — the community benchmarks
# =============================================================================
"Every column table the block A drivers dump, by run id prefix."
const A_DUMPS = Dict(
    "A1.love"  => ["love_elastic", "love_fluid"],
    "A2.disc"  => ["disc_elastic_profile", "disc_ve_profile", "disc_ve_centre_series"],
    "A3.caseA" => ["martinec_A"],
    "A5.lvz"   => ["lvz_homog_profile", "lvz_lvz_profile", "lvz_homog_series", "lvz_lvz_series"],
)

"Dump tables whose names depend on the run (the SLE cases and the rotation bodies)."
function a_dump_names(r)
    names = get(A_DUMPS, r.id, String[])
    if startswith(r.id, "A3.")
        case = replace(r.id, "A3." => "")
        append!(names, ["sle_$(case)_$(f)_profile" for f in ("fig10", "fig11", "fig12", "fig13")])
        push!(names, "sle_$(case)_history")
    end
    if startswith(r.id, "A4")
        append!(names, ["rotation_$(b)_series" for b in keys(BlockA.rotation(r).bodies)])
    end
    filter(n -> BlockA.hasdump(r, n * ".txt"), names)
end

function reduce_blockA()
    say("block A — the community benchmarks")
    runs = load_runs(joinpath(RUNS, "A"))
    isempty(runs) && (@warn "no block A output under runs/A"; return)
    @info "parsing $(length(runs)) benchmark run(s)"

    path = Published.blockA_path()
    mkpath(dirname(path)); isfile(path) && rm(path)
    ds = NCDataset(path, "c")
    try
        ds.attrib["title"] = "VILMA2 validation — block A community benchmarks, reduced"
        ds.attrib["generated_by"] = "analysis/make_reduced.jl"
        ds.attrib["generated_on"] = Dates_now()
        ds.attrib["comment"] = "One group per benchmark run. `dump_*` are the driver's own " *
            "comparison tables (model against the published reference), written verbatim; " *
            "`table_*` are the numbers parsed out of the driver's stdout. The stdout and " *
            "the SLURM epilogue themselves are not published — the verdict and the job " *
            "accounting are group attributes."

        for r in runs
            g = defGroup(ds, r.id)
            r.verdict === nothing || (g.attrib["verdict"] = r.verdict)
            for (k, v) in r.slurm
                g.attrib[k] = v
            end
            for n in a_dump_names(r)
                t = BlockA.dumptable(r, n * ".txt")
                isempty(t) && continue
                putvar(g, "dump_" * n, t;
                       attrib = Dict("long_name" => "driver dump table $n.txt"), DEFLATE...)
            end
            for (k, v) in a_tables(r)
                (v === nothing || isempty(v)) && continue
                putvar(g, "table_" * k, v;
                       attrib = Dict("long_name" => "parsed from the driver's stdout"))
            end
            @info "  $(r.id): $(verdict(r))"
        end
    finally
        close(ds)
    end
    @info @sprintf("  wrote %s (%.2f MB)", relpath(path, REPO), filesize(path) / 1e6)
end

"The numbers parsed out of one benchmark driver's stdout, keyed as Published reads them."
function a_tables(r)
    t = Dict{String,Any}()
    if r.id == "A1.love"
        d = BlockA.love(r)
        t["love_fluid"] = d.fluid; t["love_elastic"] = d.elastic
    elseif r.id == "A2.disc"
        d = BlockA.disc(r)
        t["disc_centre"] = d.centre
        t["disc_u_summary"] = d.u_summary; t["disc_n_summary"] = d.n_summary
    elseif r.id == "A3.caseA"
        d = BlockA.martinec_caseA(r)
        t["caseA_fields"] = String.(d.fields); t["caseA_peaks"] = d.peaks; t["caseA_errs"] = d.errs
    elseif startswith(r.id, "A3.")
        d = BlockA.sle(r)
        t["sle_figs"] = String.(d.figs); t["sle_table"] = d.table
        t["sle_eustatic"] = d.eustatic; t["sle_ocean"] = d.ocean; t["sle_iters"] = d.iters
    elseif startswith(r.id, "A4")
        d = BlockA.rotation(r)
        bodies = sort(collect(keys(d.bodies)))
        t["rotation_bodies"] = String.(bodies)
        for b in bodies
            t["rotation_table_" * b] = d.bodies[b]
        end
        gb = sort(collect(keys(d.G)))
        t["rotation_G_bodies"] = String.(gb); t["rotation_G"] = [d.G[b] for b in gb]
        isnan(d.kTe) || (t["rotation_kTe"] = [d.kTe])
        isnan(d.k_s) || (t["rotation_ks"] = [d.k_s])
    elseif r.id == "A5.lvz"
        d = BlockA.lvz(r)
        t["lvz_t"] = d.t; t["lvz_u"] = d.u; t["lvz_amplification"] = d.amplification
        names = sort(collect(keys(d.cases)))
        if !isempty(names)
            t["lvz_case_names"] = String.(names)
            t["lvz_case_values"] = permutedims(reduce(hcat, [collect(d.cases[n]) for n in names]))
        end
    end
    t
end

# =============================================================================
# Block B — the lateral-viscosity benchmark
# =============================================================================
# The driver's protocol files are already the reduced form (a few MB of text, in
# the format the benchmark collects), so they are copied, not converted: ℓmax 128
# for all six cases (the submission), ℓmax 64 and 256 for A and B (convergence). What is
# done here is refuse a file that is not one complete run: each case is a chain
# of 20-kyr segments that append to the same files (run_block_B.sh), so a
# duplicated header, a missing epoch or a short section is exactly what a broken
# chain would leave behind.
const BLOCKB_ROOT = "blockB"
# The normal-mode Love numbers that make test A-i checkable (analysis/NormalMode.jl),
# vendored in the model repo; published beside the runs so the check needs no model.
const BLOCKB_LOVE = joinpath(dirname(REPO), "data", "benchmarks", "love_M3-L70-V01", "mod_M3-L70-V01")

"Check one case's protocol files and copy them into `dst`; false if the run is absent."
function publish_bench_case(t, f, lmax, dst)
    dir = rundir(BLOCKB_ROOT, "$(t)-i_$(f).l$(lmax)")
    pd = joinpath(dir, Published.bench_disp_name(t, f))
    pg = joinpath(dir, Published.bench_gcm_name(t, f))
    if !(isfile(pd) && isfile(pg))
        @warn "skipping $(t)-i $(f) ℓmax $lmax (no protocol files under $dir)"
        return false
    end
    for p in (pd, pg)
        nh = count(l -> startswith(l, "# VILMA2"), eachline(p))
        nh == 1 || error("$p: $nh header blocks, expected 1 (a restart appended a second?)")
    end
    d, g = Published.read_bench_table(pd), Published.read_bench_table(pg)
    ne, np = length(Published.BENCH_EPOCHS_KYR), Published.BENCH_NPT
    if size(d, 1) < ne * np
        @warn "skipping $(t)-i $(f) ℓmax $lmax (chain incomplete: $(size(d, 1)) of $(ne * np) rows)"
        return false
    end
    size(d, 1) == ne * np ||
        error("$pd: $(size(d, 1)) rows, expected $ne epochs × $np points")
    unique(d[:, 4]) ≈ Published.BENCH_EPOCHS_KYR ||
        error("$pd: epochs $(unique(d[:, 4])) are not the protocol's")
    g[:, 1] ≈ Published.BENCH_EPOCHS_KYR ||
        error("$pg: epochs $(g[:, 1]) are not the protocol's")
    all(isfinite, d) && all(isfinite, g) || error("$(t)-i $(f) ℓmax $lmax: non-finite values")
    Published.bench_lmax(pd) == lmax || error("$pd: lmax is not $lmax")
    mkpath(dst)
    cp(pd, joinpath(dst, basename(pd)); force = true)
    cp(pg, joinpath(dst, basename(pg)); force = true)
    @info "  $(t)-i $(f) ℓmax $lmax: $ne epochs × $np points, copied"
    true
end

function reduce_blockB()
    say("block B — the lateral-viscosity benchmark")
    n = 0
    for t in Published.BENCH_TESTS, f in Published.BENCH_FORCINGS
        n += publish_bench_case(t, f, Published.BENCH_LMAX, Published.blockB_dir())
    end
    for L in Published.BENCH_CONV_LMAX, t in Published.BENCH_CONV_TESTS, f in Published.BENCH_FORCINGS
        n += publish_bench_case(t, f, L, Published.bench_case_dir(L))
    end
    isfile(BLOCKB_LOVE) || error("missing $BLOCKB_LOVE (the model repo's normal-mode table)")
    cp(BLOCKB_LOVE, Published.bench_love_path(); force = true)
    @info "  copied the M3-L70-V01 normal-mode Love numbers"
    n == 0 && @warn "no block B case was published"
    n
end

# =============================================================================
# Drive
# =============================================================================
const ALL = ["D", "V2V1", "E", "A", "B"]
want = isempty(ARGS) ? ALL : uppercase.(ARGS)
bad = setdiff(want, ALL)
isempty(bad) || error("unknown block(s): $(join(bad, ", ")). Known: $(join(ALL, ", "))")

mkpath(OUT)
"D"  in want && reduce_blockD()
"V2V1" in want && reduce_v2v1()
"E"  in want && reduce_blockE()
"A"  in want && reduce_blockA()
"B"  in want && reduce_blockB()

say("data/published")
total = 0
for (root, _, files) in walkdir(OUT), f in sort(files)
    p = joinpath(root, f)
    global total += filesize(p)
    @printf("  %-44s %8.2f MB\n", relpath(p, OUT), filesize(p) / 1e6)
end
@printf("  %-44s %8.2f MB\n", "TOTAL", total / 1e6)
