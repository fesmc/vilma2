#!/usr/bin/env julia
#
# make_vilma1_visc_ongrid.jl — resample a VILMA1-format 3-D viscosity file onto the
# radial element grid of a given densi.inp.
#
#   julia --project=. experiments/make_vilma1_visc_ongrid.jl \
#        <densi.inp> <source.nc> <out.nc> [<lid_depth_m> <lid_log10max> <log10_cap>]
#
# With the three optional arguments, VILMA2's lid rule (&vilma visc3d_lid_depth,
# visc3d_lid_log10max; model 3ba3263) is applied per element: an element lying
# wholly above lid_depth whose log10 eta exceeds lid_log10max is set to log10_cap.
# The rule is recorded in the lgvisc attribute `lid_rule`, so a pipeline can tell a
# file made with it from one made without.
#
# WHY THIS IS NECESSARY. VILMA1 maps a 3-D viscosity file onto its radial elements
# BY ARRAY INDEX, not by radius: element k takes level k of the file, and levels
# beyond its element count are never read. Measured: running
# the 164-level Bagge-derived field against the 129-element M3-L70-V01 densi.inp,
# VILMA1's own `[struct_ve3d]` echo reproduces file levels 1..129 to 2e-4 dex and
# stops at r = 6193 km -- so the top 175 km, INCLUDING THE WHOLE LITHOSPHERE,
# is never seen, and the deep-mantle profile is stretched over the whole Earth.
# The surface element came out at 10^20.95 against the file's 10^29.99.
#
# This is invisible with CLIMBER-X's stock inputs because `bagge2021.nc` is
# distributed ON the stock densi.inp's element centres -- 164 levels, 164
# elements -- so index and radius agree by construction. Any densi.inp with a
# different element count silently mis-maps the field. Two runs in this project
# hit it: the disc benchmark on the Bagge profile (129 elements) and
# D-V1-3D.l170.matched (128 elements), whose 9.2-dex lithosphere error is what
# made "matching the radial structure" appear to make agreement WORSE.
#
# THE CONTRACT, therefore: a VILMA1 3-D viscosity file must carry exactly one
# level per viscoelastic element of the densi.inp it is run with, in the same
# order. This script enforces it. It also writes the element radii into the
# output, so the file records the grid it belongs to.
#
# Element geometry, as VILMA1 builds it from densi.inp: each row after the first
# (the fluid core, which carries no viscoelastic element) is
#     r_bot  r_top  dr  rho  mu
# split into round((r_top-r_bot)/dr) elements, and the element is identified by
# its CENTRE radius. Verified against the `[struct_ve3d]` echo of a run using
# the stock file: centres 3500190.909 ... 6368500, which are exactly the radius
# coordinates of bagge2021.nc.

using NCDatasets
using Printf

length(ARGS) in (3, 6) || error("usage: make_vilma1_visc_ongrid.jl <densi.inp> <source.nc> <out.nc> " *
                                 "[<lid_depth_m> <lid_log10max> <log10_cap>]")
densi, src, out = ARGS[1:3]
lid = length(ARGS) == 6 ? parse.(Float64, ARGS[4:6]) : nothing

# ---- element centres from densi.inp -----------------------------------------
rows = Tuple{Float64,Float64,Float64}[]
for (i, ln) in enumerate(readlines(densi))
    i <= 2 && continue                      # G, then the flag line
    f = split(strip(ln))
    length(f) >= 3 || continue
    push!(rows, (parse(Float64, f[1]), parse(Float64, f[2]), parse(Float64, f[3])))
end
isempty(rows) && error("$densi: no layer rows parsed")
centres = Float64[]
bottoms = Float64[]                         # element bottom radii, for the lid rule
for (k, (rb, rt, dr)) in enumerate(rows)
    k == 1 && continue                      # fluid core: no viscoelastic element
    # VILMA1 truncates: the M3-L70 row 3480-5701 km at dr = 40 km gives 55.5 ->
    # 55 elements of 40381.8 m, matching its own [struct_ve3d] echo.
    n = max(1, floor(Int, (rt - rb) / dr))
    h = (rt - rb) / n
    append!(centres, [rb + (j - 0.5) * h for j in 1:n])
    append!(bottoms, [rb + (j - 1) * h for j in 1:n])
end
r_surface = rows[end][2]
@printf("%s: %d viscoelastic elements, centres %.1f .. %.1f km\n",
        basename(densi), length(centres), centres[1] / 1e3, centres[end] / 1e3)

# ---- source field ------------------------------------------------------------
ds = NCDataset(src)
rname = haskey(ds, "radius") ? "radius" : "r"
vname = haskey(ds, "lgvisc") ? "lgvisc" : "eta"
lon = Float64.(ds[rname == "radius" ? "lon" : "lon"][:]); lat = Float64.(ds["lat"][:])
rs  = Float64.(ds[rname][:]); a = Float64.(ds[vname][:, :, :])   # (lon, lat, r), log10 eta
close(ds)
@printf("%s: %s(%d,%d,%d), radii %.1f .. %.1f km\n", basename(src), vname,
        size(a, 1), size(a, 2), size(a, 3), rs[1] / 1e3, rs[end] / 1e3)
length(rs) == length(centres) &&
    @warn "source already has one level per element -- resampling is a no-op only if the radii agree"

"Linear in log10(eta) against radius, clamped to the source range."
function interp_col(col, rs, r)
    r <= rs[1]   && return col[1]
    r >= rs[end] && return col[end]
    k = searchsortedlast(rs, r)
    w = (r - rs[k]) / (rs[k+1] - rs[k])
    (1 - w) * col[k] + w * col[k+1]
end

nlon, nlat, nr = size(a, 1), size(a, 2), length(centres)
b = Array{Float64}(undef, nlon, nlat, nr)
for j in 1:nlat, i in 1:nlon
    col = @view a[i, j, :]
    for k in 1:nr
        b[i, j, k] = interp_col(col, rs, centres[k])
    end
end

# VILMA2's lid rule on b, in place; returns the rule as recorded in `lid_rule`.
# A function, not a top-level loop: at script top level the counter would be a
# new local of the loop (Julia soft scope) and the global never incremented.
function apply_lid!(b, bottoms, r_surface, lid_depth, lid_max, cap)
    inlid = bottoms .>= r_surface - lid_depth - 1e-6 * r_surface   # wholly above the lid base
    ncap = 0
    for k in findall(inlid), j in axes(b, 2), i in axes(b, 1)
        if b[i, j, k] > lid_max
            b[i, j, k] = cap;  ncap += 1
        end
    end
    @printf("lid rule: %d elements in the lid, %d cells capped (%.1f %% of lid cells)\n",
            count(inlid), ncap, 100 * ncap / (count(inlid) * size(b, 1) * size(b, 2)))
    return @sprintf("elements wholly above %.1f km with log10 eta > %.2f set to %.2f",
                    lid_depth / 1e3, lid_max, cap)
end

lid_rule = lid === nothing ? "none" : apply_lid!(b, bottoms, r_surface, lid...)

isfile(out) && rm(out)
d = NCDataset(out, "c")
defDim(d, "lon", nlon); defDim(d, "lat", nlat); defDim(d, "radius", nr)
defVar(d, "lon", lon, ("lon",), attrib = ["units" => "degrees_east",  "axis" => "X"])
defVar(d, "lat", lat, ("lat",), attrib = ["units" => "degrees_north", "axis" => "Y"])
defVar(d, "radius", centres, ("radius",), attrib = ["units" => "m", "axis" => "Z",
       "comment" => "element CENTRES of $(basename(densi)); VILMA1 indexes this axis, " *
                    "so it must have one level per viscoelastic element, in order"])
defVar(d, "lgvisc", b, ("lon", "lat", "radius"), attrib = [
    "long_name" => "log10 Viscosity", "units" => "Pa s",
    "comment"   => "Resampled from $(basename(src)) onto the element grid of " *
                   "$(basename(densi)) by make_vilma1_visc_ongrid.jl. VILMA1 maps this " *
                   "axis to its elements BY INDEX, not by radius.",
    "lid_rule"  => lid_rule])
close(d)
@printf("wrote %s: lgvisc(%d,%d,%d), surface level log10(eta) = %.3f\n",
        out, nlon, nlat, nr, b[1, 1, end])
