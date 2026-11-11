#!/usr/bin/env julia
#
# sle_budget.jl — decompose VILMA2's SLE grounded-ice mass source
#
#
#   julia --project=. analysis/sle_budget.jl [run/out.nc ...]
#
# With no arguments it decomposes every completed block D VILMA2 run. The
# arithmetic lives in analysis/Budget.jl, which analysis/make_reduced.jl also
# uses to build the budget series it publishes, so the table printed here and
# the curve the block D figure plots cannot drift apart.
#
# WHAT THIS MEASURES, AND WHY IT IS NOT CIRCULAR
#
# The SLE conserves whatever melt source it is handed, to machine precision, so
# its own `mass_resid` cannot tell you whether that source is right. This script
# rebuilds the source from the run's OUTPUT FIELDS (h_ice_eq, z_bed_eq, h_ice,
# C_ocean) two ways and compares both against a third quantity that the SLE
# never sees: the above-flotation volume of the ice actually lost.
#
#   as coded    -(rho_i/rho_w) * <(I - I0) * (1 - C)>
#   as required -(rho_i/rho_w) * < I*(1 - C) - I0*(1 - C0) >
#   difference  -(rho_i/rho_w) * < I0 * (C - C0) >
#
# The two agree identically when the reference is ice-free (I0 = 0) -- every
# benchmark in the suite -- and diverge on a cell that carries grounded marine
# reference ice, loses it, and becomes ocean. C0 is rebuilt here exactly as
# vilma_sle.f90 builds it, `ocean_function(topo0, ice - d_ice, C0)`, so the "as
# coded" column reproduces the run's reported sea level rather than approximating
# it: coded + flooding term == ocean-mean rsl, to ~1e-3 m.
#
# DENSITIES are the model's own (vilma_constants.f90), not oceanographic values.
# 931/1000 puts the above-flotation deglacial target at 97.3 m; 910/1028 gives
# 93.0 m, which is what sessions 26-27 quoted and is NOT the number the code
# should be measured against.
#
# QUADRATURE is the trapezoid rule in sin(lat) -- the model grid is Gaussian and
# descending in latitude, and plain cos(lat) weights are close but not exact.
# Checked against the run's own ocean fraction.
#
# SIDE OF THE PIPELINE. This is a `runs/`-side diagnostic and is deliberately
# NOT part of the figure pipeline: it prints the full budget decomposition for
# one run's out.nc, which is more than any figure needs and needs the raw file.
# The figures get the budget from data/published/, where make_reduced.jl stores
# it as a series. Run this when investigating a budget, not when drawing.

using Printf

include(joinpath(@__DIR__, "Budget.jl"))
using .Budget

const REPO = dirname(@__DIR__)

function report(path; nrow = 12)
    b = budget_series(path; nrow = nrow)
    n = length(b.time)
    println("\n", "="^76); println(replace(path, REPO * "/" => ""))
    @printf("  %d rows, %.0f -> %.0f a;  reference ocean fraction %.4f\n",
            n, b.time[1], b.time[end], b.ref_ocean_frac)
    @printf("\n%9s %9s %9s %9s %9s %9s\n",
            "time[a]", "rsl_mean", "coded", "required", "deficit", "target")
    for k in 1:n
        @printf("%9.0f %9.2f %9.2f %9.2f %9.2f %9.2f\n",
                b.time[k], b.rsl_mean[k], b.coded[k], b.required[k],
                b.coded[k] - b.required[k], b.target[k])
    end

    a = budget_at(path)
    println()
    @printf("  present day, barystatic metres (ocean fraction %.4f)\n", a.ocean_frac)
    # Whichever of the two forms the run actually used will equal the reported
    # ocean-mean rsl exactly. On a model built after the marine-reference mask
    # fix that is the `required` row; before it, the `coded` row.
    @printf("    coded    + flooding = %8.2f\n", a.coded_plus_flooding)
    @printf("    required + flooding = %8.2f   <- matches whichever form the run used\n",
            a.required_plus_flooding)
    @printf("    the run reports       %8.2f\n", a.delivered)
    @printf("    above-flotation target      = %8.2f   (delivered/target %.1f %%)\n",
            a.target, 100a.closure)
    @printf("    mask term -(ri/rw)<I0(C-C0)> = %8.2f   (what the fix recovered)\n",
            a.mask_term)
end

paths = isempty(ARGS) ?
    filter(isfile, [joinpath(REPO, "runs", "blockD_v2", d, "out.nc")
                    for d in ("D-V2-3D.l170.fe", "D-V2-1Db.l170.fe")]) : ARGS
isempty(paths) && error("no block D VILMA2 output found; pass a path explicitly")
for p in paths
    report(p)
end
