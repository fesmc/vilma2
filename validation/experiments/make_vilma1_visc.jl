#!/usr/bin/env julia
#
# make_vilma1_visc.jl — build the 3-D viscosity file VILMA1 reads, from the SAME
# field VILMA2 reads, with the SAME clamp applied.
#
#   julia --project=../analysis experiments/make_vilma1_visc.jl
#
# WHY THIS EXISTS. The two backends do not share a viscosity path. VILMA2
# reads `visc_3d_file` (model `input/bagge2021.nc`) and then applies its own
# `visc_min`/`visc_max` clamp -- reported in the run header as
# `clamp=[19.50,30.00]` -- plus an optional perturbation (`f_visc_sd`,
# `f_visc_rel`). VILMA1 reads `vilma_v1_visc_3d_file` **directly**, through its own
# `io_nc3in` unit; none of VILMA2's clamping or perturbation reaches it.
# So "point both at the Bagge field" does NOT give both backends the same field,
# and a V2-V1 difference could be a viscosity difference rather than a solver one.
#
# This script closes that gap: it applies VILMA2's clamp to VILMA2's
# own vendored input and writes the result in the layout VILMA1 expects.
#
# VERIFIED, 2026-09-19: the output agrees with CLIMBER-X's
# `visc3d_Bagge2021_min_19.5.nc` to 0.0039 dex over 12 levels spanning the
# mantle -- exactly the 0.01-dex quantization of the vendored file -- with zero
# cells differing by more than 0.01 dex. So climber's file happened to be the
# right one. That was luck, not design, and it is a 172 MB file outside this
# repository; this script makes the provenance explicit and reproducible.
#
# `f_visc_sd = 0` in every run here, which makes `f_visc_rel` inert -- it only
# scales a perturbation. If a run ever sets `f_visc_sd /= 0`, this script does
# NOT reproduce that and the comparison stops being like-for-like.
#
# The output is ~172 MB and is NOT committed; regenerate it with this script.

using NCDatasets
using Printf

const REPO  = dirname(@__DIR__)
const SRC   = get(ENV, "VISC_SRC", joinpath(REPO, "data", "input", "bagge2021.nc"))
const OUT   = get(ENV, "VISC_OUT", joinpath(REPO, "data", "reference", "bagge2021_vilma2_clamp.nc"))
# VILMA2's p%visc_min / p%visc_max, in log10(Pa s).
const FLOOR = parse(Float64, get(ENV, "VISC_FLOOR", "19.5"))
const CEIL  = parse(Float64, get(ENV, "VISC_CEIL",  "30.0"))

mkpath(dirname(OUT))
isfile(OUT) && rm(OUT)

src = NCDataset(SRC)
lon = Float64.(src["lon"][:]); lat = Float64.(src["lat"][:]); r = Float64.(src["r"][:])
@info @sprintf("source %s: %d x %d x %d", SRC, length(lon), length(lat), length(r))

d = NCDataset(OUT, "c")
defDim(d, "lon", length(lon)); defDim(d, "lat", length(lat)); defDim(d, "radius", length(r))
defVar(d, "lon",    lon, ("lon",),    attrib = ["units" => "degrees_east",  "axis" => "X"])
defVar(d, "lat",    lat, ("lat",),    attrib = ["units" => "degrees_north", "axis" => "Y"])
defVar(d, "radius", r,   ("radius",), attrib = ["units" => "km",            "axis" => "Z"])
# Name and dimension order follow CLIMBER-X's visc3d_Bagge2021.nc, which is what
# VILMA1's reader expects: lgvisc(radius, lat, lon) in the file header.
v = defVar(d, "lgvisc", Float64, ("lon", "lat", "radius"), attrib = [
    "long_name" => "log10 Viscosity", "units" => "Pa s", "grid_type" => "gaussian",
    "comment"   => "Bagge et al. (2021) as vendored in VILMA2 input/bagge2021.nc, " *
                   "clamped to [$(FLOOR),$(CEIL)] dex -- the same clamp the native solver " *
                   "applies -- so both backends integrate the same field."])

"Copy `src` eta into `v` level by level, clamped to [FLOOR, CEIL]; return the number of cells clamped."
function clamp_levels!(v, src, nr)
    nclamp = 0
    for k in 1:nr
        x = Float64.(src["eta"][:, :, k])
        nclamp += count(<(FLOOR), x) + count(>(CEIL), x)
        v[:, :, k] = clamp.(x, FLOOR, CEIL)
    end
    return nclamp
end

nclamp = clamp_levels!(v, src, length(r))
close(d); close(src)
@info @sprintf("wrote %s (%d cells clamped)", OUT, nclamp)
