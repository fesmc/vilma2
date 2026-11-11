#!/usr/bin/env julia
#
# make_vilma1_densi.jl — a densi.inp for VILMA1 built from VILMA2's own
# radial layer table, so the two backends share a radial structure.
#
#   julia --project=. experiments/make_vilma1_densi.jl [out.inp]
#
# WHY. The V2-V1 comparison matched the viscosity field (make_vilma1_visc.jl) and
# showed VILMA1's time step to be converged, yet a difference remains. The one
# input still unmatched is the RADIAL STRUCTURE, and it is unmatched in a way
# that plausibly explains the signature we see -- VILMA1 lagging on the FAST
# response and agreeing on the slow one, which is what a different lithosphere
# or shallow discretisation does.
#
# WHAT ACTUALLY DIFFERS. The element-size schedules already agree: both use the
# VEGA 5/10/40 km spacing (VILMA2's DR_LITHO/DR_UPPER/DR_LOWER, with
# breakpoints at 70 and 670 km depth), and both end up near 220 radial elements.
# What differs is the LAYER TABLE:
#
#   VILMA2   11 PREM layers; interfaces at 3480, 3630, 4287, 4943, 5600,
#                 5701, 5771, 5971, 6151, 6291 km. Per-layer PREM rho and mu.
#                 Elastic lithosphere 6291-6371 km, i.e. 80 km.
#   VILMA1 stock   7 layers; interfaces at 3480, 5701, 5951, 6271, 6311, 6366 km.
#                 Constant rho = 4449.4 and mu = 1.4519e11 / 0.67e11 through the
#                 mantle. Lithosphere boundary at 6271 km, i.e. 100 km.
#
# Only ONE interface (5701 km) is common, and the lithospheric boundary differs
# by 20 km -- which is exactly the parameter the short-wavelength, fast response
# is most sensitive to.
#
# The stock rho/mu columns look like placeholders, consistent with vg%l_prem = 1
# overriding them from a polynomial PREM. So this file is only fully effective
# with **vilma_v1_l_prem = 0**, which is why that was made a namelist parameter.
# Running both l_prem settings against the same file also MEASURES whether the
# rho/mu columns are used at all -- if the two runs agree, they are not.
#
# FORMAT (list-directed; radii in metres, ascending):
#   line 1   gravitational constant G
#   line 2   five integer flags, copied from the stock file unchanged
#   line 3+  r_bot  r_top  dr  rho  mu
#
# THE CORE. VILMA2 gives the core mu = 0 (fluid, rheology = 2); the stock
# densi.inp carries 5.4519e10. Written as 0 here, which is the physically right
# value and what VILMA2 integrates. If VILMA1 rejects it, set
# CORE_MU=5.4519e10 in the environment to fall back to the stock value -- and
# then say so, because the two backends no longer share a core.

using Printf

const REPO  = dirname(@__DIR__)
const STOCK = joinpath(get(ENV, "VILMA1_INPUT", joinpath(REPO, "experiments", "vilma1_input")), "densi.inp")
const OUT   = length(ARGS) >= 1 ? ARGS[1] :
              joinpath(REPO, "experiments", "vilma1_input_vilma2grid", "densi.inp")

# VILMA2's 11-layer PREM table, exactly as the block D namelist supplies it
# (experiments/timing/blockD.nml.tmpl), surface-first.
const R_TOP = [6371.0, 6291.0, 6151.0, 5971.0, 5771.0, 5701.0, 5600.0, 4943.0, 4287.0, 3630.0, 3480.0] .* 1e3
const R_BOT = [6291.0, 6151.0, 5971.0, 5771.0, 5701.0, 5600.0, 4943.0, 4287.0, 3630.0, 3480.0,    0.0] .* 1e3
const RHO   = [3175.1, 3367.3, 3489.4, 3848.4, 3984.9, 4411.9, 4628.0, 4986.1, 5318.7, 5538.8, 10987.2]
const MU    = [5.727e10, 6.655e10, 7.732e10, 1.051e11, 1.225e11, 1.638e11, 1.926e11,
               2.313e11, 2.690e11, 2.913e11, 0.0]

const R_EARTH     = 6371.0e3
const DR_LITHO    = 5.0e3            # VILMA2 vilma_radial_fe.f90
const DR_UPPER    = 10.0e3
const DR_LOWER    = 40.0e3
const DEPTH_LITHO = 70.0e3
const DEPTH_UPPER = 670.0e3
const CORE_MU     = parse(Float64, get(ENV, "CORE_MU", "0.0"))

"VILMA2's dr_target: element size from mid-layer depth."
function dr_target(depth)
    depth <= DEPTH_LITHO ? DR_LITHO : depth <= DEPTH_UPPER ? DR_UPPER : DR_LOWER
end

# Preserve the stock header: G on line 1, the flag line on line 2. Neither is
# ours to invent -- G is a physical constant VILMA1 reads, and the flags select
# code paths we have no documentation for.
stock = readlines(STOCK)
length(stock) >= 2 || error("$STOCK too short to carry a header")
gconst, flags = strip(stock[1]), strip(stock[2])
@info "header taken from the stock file: G = $gconst, flags = $flags"

open(OUT, "w") do io
    println(io, gconst)
    println(io, flags)
    ne = 0
    # VILMA1 lists layers innermost-first (ascending r), VILMA2 surface-first.
    for i in length(R_TOP):-1:1
        r0, r1 = R_BOT[i], R_TOP[i]
        dr = dr_target(R_EARTH - 0.5*(r0 + r1))
        mu = (MU[i] == 0.0) ? CORE_MU : MU[i]
        n  = max(1, ceil(Int, (r1 - r0)/dr))
        ne += n
        @printf(io, "%12.4f  %12.4f  %10.4f  %12.4f  %14.6e\n", r0, r1, dr, RHO[i], mu)
    end
    @printf("wrote %s: %d layers, %d elements (VILMA2 reports 220)\n",
            OUT, length(R_TOP), ne)
end
