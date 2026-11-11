#!/usr/bin/env julia
#
# make_uniform_visc.jl — a LATERALLY UNIFORM viscosity field, in both backends'
# file formats, from one set of numbers.
#
#   julia --project=../analysis experiments/make_uniform_visc.jl
#
# WHY. VILMA2 and VILMA1 disagree on the 3-D problem by 12.6 m rms in rsl
# (21 ka) and 2.70 mm/yr in present-day uplift, from a verified-identical initial
# state, with the viscosity field made identical and VILMA1's time step shown to
# be converged (|V(24) - V(48)| = 0.006 m rms against a 12.6 m difference). Two
# candidates remain:
#
#   (a) the RADIAL structures are not matched. VILMA2 uses the explicit
#       11-layer PREM from its namelist on 220 radial elements with a 1e40 Pa s
#       lithosphere; VILMA1 builds its own from polynomial PREM plus densi.inp,
#       taking its lithosphere from the 3-D file. Same nominal Earth, different
#       discretisation.
#   (b) a genuine difference in how lateral viscosity is treated.
#
# Running both against a laterally UNIFORM field separates them. With no lateral
# structure there is nothing for (b) to act on, so any residual difference is
# (a) -- or a difference in the 1-D response itself, which the community
# benchmarks in block A already constrain for VILMA2. If they agree here and
# disagree with the Bagge field, the difference is specifically in the lateral
# treatment, which is the part with no community benchmark at all.
#
# This is VILMA2's own A6 consistency check ("uniform eta supplied as a 3-D
# field must reproduce the 1-D path") extended across codes.
#
# THE MEAN. Per radial level, the unweighted arithmetic mean of log10(eta) --
# i.e. the geometric mean of eta. That is exactly what VILMA2 does when it
# demotes an element (vilma_coupling.f90: visc_unif(:,r) = sum(visc_node(:,r))/nh),
# and what visc3d_tol = 99 produces for structure 1Db. Area weighting is
# deliberately NOT applied, to match that convention. It does not affect the
# test either way: both backends read the SAME numbers, so the comparison is
# clean whatever the mean is.
#
# Two files, one source of truth:
#   data/reference/bagge2021_uniform_vilma2.nc     eta(r,lat,lon),     for visc_3d_file
#   data/reference/bagge2021_uniform_vilma1.nc  lgvisc(radius,...), for vilma_v1_visc_3d_file

using NCDatasets
using Printf
using Statistics

const REPO  = dirname(@__DIR__)
const SRC   = joinpath(REPO, "data", "reference", "bagge2021_vilma2_clamp.nc")
const OUT_F = joinpath(REPO, "data", "reference", "bagge2021_uniform_vilma2.nc")
const OUT_V = joinpath(REPO, "data", "reference", "bagge2021_uniform_vilma1.nc")

isfile(SRC) || error("$SRC missing -- run experiments/make_vilma1_visc.jl first")
src = NCDataset(SRC)
lon = Float64.(src["lon"][:]); lat = Float64.(src["lat"][:]); r = Float64.(src["radius"][:])
nlon, nlat, nr = length(lon), length(lat), length(r)

prof = Vector{Float64}(undef, nr)
spread = Vector{Float64}(undef, nr)
for k in 1:nr
    x = Float64.(src["lgvisc"][:, :, k])
    prof[k]   = mean(x)                       # geometric mean of eta
    spread[k] = maximum(x) - minimum(x)       # lateral range being removed [dex]
end
close(src)
@printf("radial profile: log10(eta) in [%.3f, %.3f]; lateral range removed up to %.3f dex\n",
        minimum(prof), maximum(prof), maximum(spread))

function write_field(path, varname, dimnames, radname, T)
    mkpath(dirname(path)); isfile(path) && rm(path)
    d = NCDataset(path, "c")
    defDim(d, "lon", nlon); defDim(d, "lat", nlat); defDim(d, radname, nr)
    defVar(d, "lon", lon, ("lon",), attrib = ["units" => "degrees_east",  "axis" => "X"])
    defVar(d, "lat", lat, ("lat",), attrib = ["units" => "degrees_north", "axis" => "Y"])
    defVar(d, radname, r, (radname,), attrib = ["units" => "km", "axis" => "Z"])
    v = defVar(d, varname, T, dimnames, attrib = [
        "long_name" => "log10 Viscosity (laterally uniform)", "units" => "Pa s",
        "comment"   => "Per-level geometric mean of the clamped Bagge (2021) field. " *
                       "Laterally uniform BY CONSTRUCTION, so VILMA2 and VILMA1 " *
                       "integrate the same radial profile and any difference is radial " *
                       "structure, not lateral treatment."])
    for k in 1:nr
        v[:, :, k] = fill(T(prof[k]), nlon, nlat)
    end
    close(d)
    @info "wrote $path"
end

# VILMA2 reads `eta(r, lat, lon)` as float32 (model input/bagge2021.nc layout).
write_field(OUT_F, "eta",    ("lon", "lat", "r"),      "r",      Float32)
# VILMA1 reads `lgvisc(radius, lat, lon)` as float64 (CLIMBER-X layout).
write_field(OUT_V, "lgvisc", ("lon", "lat", "radius"), "radius", Float64)
