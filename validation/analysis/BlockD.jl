"""
    BlockD

Loading helpers for the last-deglaciation maps.

The driver writes `h_ice`, `rsl`, `z_bed`, `C_ocean` and `bsl` on its own Gauss
grid, north-to-south in latitude and 0…360 in longitude. Everything here returns
fields ready to plot as **(lon, lat)**: longitude rolled to −180…180 so the
Atlantic is not split down the middle of the frame, and latitude flipped
south-to-north so Makie's y-axis runs the way a map should.

Sign conventions, because they are easy to get backwards:

  * `rsl` is the FULL relative-sea-level change field (N − u + Δφ) against the
    run's reference state, defined everywhere, not just over ocean. The paper's
    "anomaly relative to present day" is `rsl(t) − rsl(0)`, which is what
    `rsl_anomaly` returns — the model's own present-day state is the datum, not
    the reference topography, so a non-zero model PD bias cancels.
  * `z_bed = z_bed_eq − rsl`, so bedrock UPLIFT rate is `−d(rsl)/dt`.
    `uplift_rate` returns mm/yr, positive up.

SIDE OF THE PIPELINE. This module reads `runs/` and is therefore part of the
REDUCTION half: it is included by `analysis/make_reduced.jl` and by nothing
else. The figures read `data/published/` through `analysis/Published.jl`, which
serves the same interface from the reduced files. If you are writing a figure
and reach for this module, you want Published instead.
"""
module BlockD

using NCDatasets
using Statistics

export DeglacRun, load_deglac, slice_at, rsl_anomaly, uplift_rate, srf_elevation,
       nearest_t, ocean_mean, polar_motion, regrid, site_series, lat_weights

struct DeglacRun
    id::String
    dir::String
    lmax::Int
    visc::String          # "1d" | "3d"
    lon::Vector{Float64}  # −180…180, ascending
    lat::Vector{Float64}  # south→north, ascending
    time::Vector{Float64} # yr, negative = past
    ds::NCDataset
    roll::Int             # longitude shift applied
end

Base.show(io::IO, r::DeglacRun) =
    print(io, "DeglacRun($(r.id), ℓmax $(r.lmax), $(r.visc), $(length(r.time)) steps)")

"""
    load_deglac(dir; id) -> DeglacRun

Open `<dir>/out.nc` lazily (slices are read on demand — the full files are ~500 MB
each and only a handful of time levels are ever plotted).
"""
function load_deglac(dir::AbstractString; id::AbstractString = basename(dir))
    path = joinpath(dir, "out.nc")
    isfile(path) || return nothing
    ds = NCDataset(path)
    lon = Float64.(ds["lon"][:])
    lat = Float64.(ds["lat"][:])
    t   = Float64.(ds["time"][:])

    # 0…360 → −180…180. circshift by the number of cells at or above 180°.
    roll = count(>=(180.0), lon)
    lonr = circshift(lon, roll)
    lonr[lonr .>= 180.0] .-= 360.0

    m = match(r"\.?l(\d+)", id)
    lmax = m === nothing ? 0 : parse(Int, m[1])
    # Run ids spell the structure `3D` and `1Db`, so the test has to be
    # case-insensitive AND has to reject `1Db` before it matches `D`. Getting
    # this wrong is silent: every run reports itself as 1-D and the figures
    # still draw.
    visc = occursin(r"1D"i, id) ? "1d" : occursin(r"3D"i, id) ? "3d" : "1d"
    DeglacRun(id, dir, lmax, visc, lonr, reverse(lat), t, ds, roll)
end

"Index of the time level nearest `yr`."
nearest_t(r::DeglacRun, yr) = argmin(abs.(r.time .- yr))

"""
Apply the run's longitude roll and latitude flip to one slab.

NCDatasets reports dimensions in the reverse of the netCDF header order, so a
variable declared `(time, lat, lon)` in the file is `(lon, lat, time)` here and a
time slice arrives already as **(lon, lat)** — the layout Makie's
`heatmap!(ax, x, y, z)` wants. No transpose is needed; only the roll (dim 1) and
the south-to-north flip (dim 2).
"""
function orient(r::DeglacRun, a::AbstractMatrix)
    b = circshift(a, (r.roll, 0))     # (lon, lat): roll longitude
    reverse(b, dims = 2)              # lat came north→south; flip to south→north
end

"""
    slice_at(r, name, yr) -> Matrix (lon, lat)

One time level of a 3-D variable, nearest to `yr`.
"""
slice_at(r::DeglacRun, name::AbstractString, yr::Real) =
    orient(r, Float64.(r.ds[name][:, :, nearest_t(r, yr)]))

"""
    rsl_anomaly(r, yr) -> Matrix [m]

Relative sea level at `yr` minus the model's own present-day field. Positive =
sea level higher then than now, which over a forebulge reads as land that has
since risen.
"""
rsl_anomaly(r::DeglacRun, yr::Real) = slice_at(r, "rsl", yr) .- slice_at(r, "rsl", 0.0)

"""
    uplift_rate(r; yr = 0.0, window = 500.0) -> Matrix [mm/yr]

Present-day bedrock uplift rate, centred differences over `window` years where
both endpoints exist. z_bed = z_bed_eq − rsl, so this is −Δrsl/Δt. A 500 yr
window (five coupling steps) rather than one step: the SLE is converged only to
its own tolerance each step, and differencing adjacent steps amplifies that into
visible speckle without changing the signal.
"""
function uplift_rate(r::DeglacRun; yr::Real = 0.0, window::Real = 500.0)
    i1 = nearest_t(r, yr)
    i0 = nearest_t(r, yr - window)
    dt = r.time[i1] - r.time[i0]
    dt == 0 && error("uplift_rate: zero time window at yr=$yr")
    a1 = orient(r, Float64.(r.ds["rsl"][:, :, i1]))
    a0 = orient(r, Float64.(r.ds["rsl"][:, :, i0]))
    @. -(a1 - a0) / dt * 1000.0          # m/yr → mm/yr, positive up
end

"""
    srf_elevation(r, yr; h_min = 10.0) -> (z_srf, ice_only)

Surface elevation z_bed + h_ice, and the same field masked to NaN wherever the
ice is thinner than `h_min`. The masked copy is what gets its own colour scale
and contours on top of the sea-level field; NaN renders transparent in Makie.
"""
function srf_elevation(r::DeglacRun, yr::Real; h_min::Real = 10.0)
    h = slice_at(r, "h_ice", yr)
    z = slice_at(r, "z_bed", yr)
    zs = z .+ h
    ice = copy(zs)
    ice[h .< h_min] .= NaN
    (zs, ice)
end

# ---------------------------------------------------------------------------
# Global scalars
# ---------------------------------------------------------------------------
"""
    lat_weights(lat) -> Vector

Quadrature weights in sin(lat), normalised to 1, for either latitude ordering.
The model grid is Gaussian, so plain cos(lat) weights are close but not exact —
close enough to hide a real error in a budget that is quoted to a tenth of a
metre, which is why this is the trapezoid rule in sin(lat) instead.
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

`<C·rsl> / <C>` at the time level nearest `yr`: relative sea level averaged over
the ocean of that moment, which is the barystatic quantity the sea-level budget
is judged on. Note this is NOT the `bsl` diagnostic — `bsl` is the model's own
eustatic bookkeeping, while this is an area average of the solved field.
"""
function ocean_mean(r::DeglacRun, yr::Real)
    it = nearest_t(r, yr)
    C   = Float64.(r.ds["C_ocean"][:, :, it])
    rsl = Float64.(r.ds["rsl"][:, :, it])
    w   = lat_weights(Float64.(r.ds["lat"][:]))
    nlon = length(r.lon)
    sph(F) = sum((sum(@view F[:, j]) / nlon) * w[j] for j in eachindex(w))
    sph(C .* rsl) / sph(C)
end

"""
    polar_motion(r) -> (m1, m2) [as written by the driver]

The rotation vector perturbation, `rot_m_re` + i·`rot_m_im`. Present only in
runs made on VILMA2 `b4b19d7` or later; `nothing` otherwise, so a caller
can degrade rather than crash on an older run.
"""
function polar_motion(r::DeglacRun)
    (haskey(r.ds, "rot_m_re") && haskey(r.ds, "rot_m_im")) || return nothing
    (Float64.(r.ds["rot_m_re"][:]), Float64.(r.ds["rot_m_im"][:]))
end

"""
    site_series(r, name, lat, lon) -> (values, dlat, dlon)

One grid cell's full time series, with the offset of that cell from the
requested point. Indices are resolved on the ORIENTED axes and then mapped back
to the file's own layout — the roll is circular and the latitude flip is a
reversal — because `r.ds[...]` is the raw file and knows nothing about either.
"""
function site_series(r::DeglacRun, name::AbstractString, lat::Real, lon::Real)
    i = argmin(abs.(r.lon .- lon))
    j = argmin(abs.(r.lat .- lat))
    ifile = mod1(i - r.roll, length(r.lon))
    jfile = length(r.lat) - j + 1
    (Float64.(r.ds[name][ifile, jfile, :]), abs(r.lat[j] - lat), abs(r.lon[i] - lon))
end

# ---------------------------------------------------------------------------
# Comparing runs that are not on the same grid
# ---------------------------------------------------------------------------
"""
    regrid(field, lon_src, lat_src, lon_dst, lat_dst) -> Matrix

Bilinear interpolation of an oriented (lon, lat) field onto another oriented
grid. Longitude wraps; latitude clamps at the poles.

This exists for block E1, where the truncations 64/128/256 land on four
different Gauss grids and the convergence question — does the answer stop moving
as ℓmax rises — cannot be asked cell by cell until they are on one grid. The
COARSEST grid of the pair is always the target: interpolating the coarse field
up would invent structure it does not have and report the invention as
convergence error.
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

"""
    regrid(field, src::DeglacRun, dst::DeglacRun)

`field` is on `src`'s oriented grid; return it on `dst`'s.
"""
regrid(field::AbstractMatrix, src::DeglacRun, dst::DeglacRun) =
    size(field) == (length(dst.lon), length(dst.lat)) ? field :
    regrid(field, src.lon, src.lat, dst.lon, dst.lat)

end # module
