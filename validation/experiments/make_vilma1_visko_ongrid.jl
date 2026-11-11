#!/usr/bin/env julia
#
# make_vilma1_visko_ongrid.jl — a 1-D visko.inp from a laterally uniform VILMA1
# viscosity file, on the element grid of a given densi.inp.
#
#   julia --project=../analysis experiments/make_vilma1_visko_ongrid.jl \
#        <densi.inp> <uniform.nc> <out.inp>
#
# WHY. The 1-D control of the V2-V1 comparison has to be the radial mean OF THE
# FIELD THE 3-D RUN INTEGRATES, or it is not a control: `D-V1-1Db` read the stock
# `visko.inp`, whose lower mantle is 10^21.1-10^21.6 against the 10^21.5-10^23.1
# that VILMA2's 1Db actually integrates, and whose lithosphere is 120.7 km
# against 80 km. That is a different Earth, not a control.
#
# There are two ways to hand VILMA1 a 1-D structure -- `l_visc_3d = .false.` with
# a visko.inp, and `l_visc_3d = .true.` with a laterally uniform 3-D file -- and
# the paper derives its 1-D structure from Bagge, so both should exist and should
# agree. This script builds the first FROM the second, so they are the same
# numbers by construction and any difference between the two runs is VILMA1's file
# handling rather than the Earth model.
#
# FORMAT (from the stock file and experiments/make_vilma1_earth_m3l70.jl):
# `r_bot r_top eta eta`, ascending in r, the two eta columns equal, and a first
# row for the fluid core -- which carries mu = 0 and is inviscid regardless, but
# VILMA1 wants a finite number there. One row PER ELEMENT here rather than per
# layer: the source file has one level per element already, so a per-element
# table transfers it without a resampling step to get wrong.
#
# Element geometry is VILMA1's, as make_vilma1_visc_ongrid.jl documents it: each
# densi.inp row after the first is split into floor((r_top-r_bot)/dr) elements.

using NCDatasets
using Printf
using Statistics

length(ARGS) == 3 || error("usage: make_vilma1_visko_ongrid.jl <densi.inp> <uniform.nc> <out.inp>")
densi, src, out = ARGS

rows = Tuple{Float64,Float64,Float64}[]
for (i, ln) in enumerate(readlines(densi))
    i <= 2 && continue                      # G, then the flag line
    f = split(strip(ln))
    length(f) >= 3 || continue
    push!(rows, (parse(Float64, f[1]), parse(Float64, f[2]), parse(Float64, f[3])))
end
isempty(rows) && error("$densi: no layer rows parsed")

core = rows[1]                              # fluid core: no viscoelastic element
edges = Tuple{Float64,Float64}[]
for (k, (rb, rt, dr)) in enumerate(rows)
    k == 1 && continue
    n = max(1, floor(Int, (rt - rb) / dr))
    h = (rt - rb) / n
    for j in 0:(n - 1)
        push!(edges, (rb + j*h, rb + (j + 1)*h))
    end
end

d = NCDataset(src)
vname = haskey(d, "lgvisc") ? "lgvisc" : "eta"
lg = Float64.(d[vname][:, :, :])            # (lon, lat, level)
nlev = size(lg, 3)
close(d)

nlev == length(edges) ||
    error("$src has $nlev levels but $densi has $(length(edges)) viscoelastic elements")

# The whole point is that this file is laterally uniform. Refuse anything else
# rather than silently collapsing a 3-D field to its mean by a convention this
# script would then own.
spread = maximum(maximum(lg[:, :, k]) - minimum(lg[:, :, k]) for k in 1:nlev)
spread < 1e-9 || error(@sprintf("%s is NOT laterally uniform: max spread %.3e dex", src, spread))

prof = [10.0^mean(lg[:, :, k]) for k in 1:nlev]

open(out, "w") do io
    @printf(io, "%14.6e  %14.6e  %14.6e  %14.6e\n", core[1], core[2], 1.0e19, 1.0e19)
    for (k, (rb, rt)) in enumerate(edges)
        @printf(io, "%14.6e  %14.6e  %14.6e  %14.6e\n", rb, rt, prof[k], prof[k])
    end
end
@printf("wrote %s: %d element rows + core, log10(eta) in [%.3f, %.3f], surface %.3f\n",
        out, length(edges), log10(minimum(prof)), log10(maximum(prof)), log10(prof[end]))
