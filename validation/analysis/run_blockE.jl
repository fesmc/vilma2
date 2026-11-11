#!/usr/bin/env julia
#
# run_blockE.jl — the two sensitivity axes the pruned list keeps.
#
#   julia --project=../analysis analysis/run_blockE.jl
#
#   E1  lmax 64 / 128 / 256, against the production 170 — the numerical
#       convergence of the headline result, which a reviewer will ask for.
#   E3  visc3d_tol 1e-3 / 1e-2 / 1e-1 and tol = 0 (every Maxwell element treated
#       as laterally 3-D) — VILMA2's central speed/accuracy approximation,
#       which has to be justified rather than stated.
#
# E2 (forcing cadence) is not here and is not a missing script: the coupling
# interval is the forcing file's own time axis, not a namelist parameter, so
# testing 500 yr needs a model or data change rather than a run. See
# experiments.md §E.
#
# HOW CONVERGENCE IS MEASURED, AND WHY IT IS NOT SELF-REFERENTIAL. Each axis is
# measured against its OWN most-resolved member — lmax 256 for E1, tol = 0 for
# E3 — not against the production setting. Measuring against lmax 170 would
# report the production run as exact by construction and say nothing about
# whether 170 is enough.
#
# The lmax runs land on four different Gauss grids, so the fields are brought
# onto the COARSER grid of each pair before differencing. Interpolating the
# coarse field up instead would invent structure it does not have and then
# report the invention as convergence error.
#
# THAT DIFFERENCING HAPPENS AT REDUCTION, NOT HERE. Every quantity this script
# plots is one number per run, and the numbers that are differences need the
# full fields of two runs on two different grids — which is exactly what the
# published epoch subset does not carry. So `analysis/make_reduced.jl` computes
# them once and `data/published/blockE.nc` is a table of scalars. This script
# draws that table; it opens no model output at all.

using CairoMakie
using Printf
using Statistics

include(joinpath(@__DIR__, "Published.jl"))
include(joinpath(@__DIR__, "FigStyle.jl"))
using .Published
using .FigStyle

const REPO = dirname(@__DIR__)

# E1 is the resolution ladder, E3 the lateral-promotion threshold; which runs
# feed each, and which member is the reference, is declared in
# analysis/make_reduced.jl, because that is where the comparison is made. Here
# they arrive already reduced, in plotting order.
e1, e3, meta = load_blockE()

# =============================================================================
# E1 — spectral convergence
# =============================================================================
function fig_lmax()
    length(e1) < 2 && (@warn "E1 needs at least two truncations"; return nothing)
    ref_lmax = meta.ref_lmax
    @info "E1: reference is ℓmax $ref_lmax"

    L    = [Float64(x.lmax) for x in e1]
    om   = [x.ocean_mean for x in e1]
    # The above-flotation target is rebuilt from each run's OWN reference fields,
    # so it moves when the reference is remapped onto a different Gauss grid.
    # Plotting it beside the delivered sea level is what separates "the solver
    # converges" from "the reference state shifted under it" -- the two are
    # indistinguishable on the delivered curve alone.
    tgt  = [x.target for x in e1]
    clo  = [x.closure for x in e1]
    drsl = [x.rms_rsl for x in e1]
    dupl = [x.rms_upl for x in e1]
    wall = [x.wall_min for x in e1]

    fig = Figure(size = (1250, 400))
    ax1 = Axis(fig[1, 1], xscale = log2, xticks = (L, string.(Int.(L))),
               xlabel = "ℓmax", ylabel = "sea level at present day [m]")
    lines!(ax1, L, om; color = C_3D, label = "delivered (ocean-mean rsl)")
    scatter!(ax1, L, om; color = C_3D)
    lines!(ax1, L, tgt; color = C_1D, linestyle = :dash, label = "above-flotation volume lost")
    scatter!(ax1, L, tgt; color = C_1D)
    axislegend(ax1, position = :rc)

    # Explicit decade-and-subdivision ticks: Makie's log default labels these as
    # 10^-0.5, which is not a number a reader converts in their head.
    logticks(vals) = (Float64.(vals), [v >= 1 ? string(v) : string(v) for v in vals])
    ax2 = Axis(fig[1, 2], xscale = log2, yscale = log10, xticks = (L, string.(Int.(L))),
               yticks = logticks([0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5]),
               xlabel = "ℓmax", ylabel = "rms difference from ℓmax $ref_lmax")
    keep = .!isnan.(drsl)
    lines!(ax2, L[keep], drsl[keep]; color = C_3D, label = "rsl at 21 ka [m]")
    scatter!(ax2, L[keep], drsl[keep]; color = C_3D)
    lines!(ax2, L[keep], dupl[keep]; color = C_1D, label = "uplift rate [mm/yr]")
    scatter!(ax2, L[keep], dupl[keep]; color = C_1D)
    axislegend(ax2, position = :lb)

    ax3 = Axis(fig[1, 3], xscale = log2, yscale = log10, xticks = (L, string.(Int.(L))),
               yticks = logticks([2, 5, 10, 20, 50, 100]),
               xlabel = "ℓmax", ylabel = "wall clock, full deglaciation [min]")
    kw = .!isnan.(wall)
    lines!(ax3, L[kw], wall[kw]; color = C_3D)
    scatter!(ax3, L[kw], wall[kw]; color = C_3D)

    panels!(fig, [(1, 1), (1, 2), (1, 3)])
    colgap!(fig.layout, 18)

    savefig(fig, "lmax-convergence";
        title = "Spectral convergence of the deglaciation",
        caption = """
        Convergence of the last-deglaciation solution with spectral truncation, on the
        3-D viscosity field with everything else held at the block D configuration.
        (a) Present-day sea level against ℓmax: the delivered ocean-mean relative sea
        level, and the above-flotation volume of reference ice lost, both in barystatic
        metres and both rebuilt from each run's own fields. (b) Root-mean-square
        difference of each truncation from the ℓmax $ref_lmax solution, for relative sea
        level at 21 ka and for the present-day uplift rate; fields are compared on the
        coarser grid of each pair. (c) Measured wall clock for the complete 26 kyr run at
        each truncation, on 64 threads of one Levante node.
        """,
        notes = """
        The reference is the most-resolved member, ℓmax $ref_lmax, not the production
        truncation: measuring against ℓmax 170 would make the production run exact by
        construction and answer nothing. The production choice of 170 is what (b) and (c)
        together are meant to justify — how much accuracy is left on the table, against
        what it would cost to collect it.
        Differencing across truncations requires a common grid, so the finer field is
        interpolated bilinearly down onto the coarser one. Interpolating the coarse field
        up would invent structure it does not carry and report that as convergence error.
        The ℓmax 64/128/256 runs are `runs/blockE`; the 170 point is the block D
        production run itself, so (a) and (b) place the headline result on this curve
        rather than a re-run that resembles it.
        THE ℓmax 170 POINT SITS ABOUT A METRE BELOW THE OTHER THREE IN (a), AND IS
        CORRESPONDINGLY HIGH IN (b). That is not a truncation error. The LGM reference
        state is stored on the ℓmax 170 Gauss grid, so the 170 run reads it natively while
        every other truncation interpolates it; interpolation moves the above-flotation
        volume of the reference ice, and both curves in (a) move together by the same
        amount. The delivered/target ratio is $(join([@sprintf("%.2f", c) for c in clo], " / ")) %
        at ℓmax $(join(Int.(L), " / ")) — flat to about a tenth of a percent — which is the
        invariant that says the solver is converging and the datum is what moved. Read the
        clean convergence statement off the 64 → 128 → 256 sequence, which shares one
        interpolation path.
        """,
        source = @__FILE__,
        runs = [Published.blockE_path()])

    for (i, x) in enumerate(e1)
        record!("blockE", "E1_l$(x.lmax)_ocean_mean_m" => round(om[i], digits = 2))
        isnan(drsl[i]) || record!("blockE",
            "E1_l$(x.lmax)_rms_vs_ref_m" => round(drsl[i], digits = 2),
            "E1_l$(x.lmax)_uplift_rms_vs_ref_mmyr" => round(dupl[i], digits = 3))
    end
    record!("blockE", "E1_reference_lmax" => ref_lmax)
    for (i, x) in enumerate(e1)
        record!("blockE", "E1_l$(x.lmax)_target_m"      => round(tgt[i], digits = 2),
                          "E1_l$(x.lmax)_closure_pct"   => round(clo[i], digits = 2))
    end
    e1
end

# =============================================================================
# E3 — visc3d_tol, the lateral-promotion threshold
# =============================================================================
function fig_tol()
    length(e3) < 2 && (@warn "E3 needs at least two tolerances"; return nothing)
    @info "E3: reference is $(meta.ref_tol_id)"

    x  = 1:length(e3)
    n3 = [v.ne3d for v in e3]
    nt = meta.ne_tot
    drsl = [v.rms_rsl for v in e3]
    dupl = [v.rms_upl for v in e3]
    wall = [v.wall_min for v in e3]

    # Phrase the accuracy result from the MEASURED numbers rather than asserting it,
    # so that if a future change makes the approximation cost something the caption
    # says so instead of repeating that it is free.
    fin   = [d for d in drsl if !isnan(d)]
    finu  = [d for d in dupl if !isnan(d)]
    mx    = isempty(fin)  ? NaN : maximum(fin)
    mxu   = isempty(finu) ? NaN : maximum(finu)
    nzero = count(d -> !isnan(d) && d == 0.0, drsl)
    mx2   = isempty(fin) ? NaN : minimum(fin)     # the tight-threshold floor
    exact = nzero > 0
    acc_phrase = exact ?
        "$(nzero) of the $(length(fin)) thresholds are bit-identical to it, and the " *
        "largest departure anywhere in the range is " * @sprintf("%.3f", mx) *
        " m rms in relative sea level at 21 ka and " * @sprintf("%.4f", mxu) *
        " mm/yr in present-day uplift." :
        "the largest departure over the whole range is " * @sprintf("%.3f", mx) *
        " m rms in relative sea level at 21 ka and " * @sprintf("%.4f", mxu) *
        " mm/yr in present-day uplift."
    acc_short = exact ?
        "an answer that is bit-identical over two decades of the threshold and departs " *
        "by at most " * @sprintf("%.3f", mx) * " m rms across the whole range." :
        @sprintf("%.3f", mx) * " m rms in relative sea level at 21 ka."

    fig = Figure(size = (860, 400))
    labs = [v.label for v in e3]

    ax1 = Axis(fig[1, 1], xticks = (collect(x), labs),
               xlabel = "visc3d_tol [dex of lateral spread]",
               ylabel = "radial elements advanced as 3-D")
    barplot!(ax1, collect(x), Float64.(n3); color = C_3D, width = 0.5)
    nt > 0 && hlines!(ax1, [Float64(nt)], color = C_GREY, linestyle = :dash)

    # NO ACCURACY PANEL. It plotted the rms difference from tol = 0 on a log axis,
    # and those differences are now exactly zero -- 1e-3 and 1e-2 are bit-identical
    # to tol = 0 -- so the panel rendered empty. Before the degree-1 Z6 defect was
    # fixed it showed 0.403 m at every threshold, which read as the price of the
    # approximation and was nothing of the kind. A panel of
    # zeros invites the reader to hunt for structure that is not there; the result
    # is stated in the caption instead, from the measured numbers, and both are
    # still recorded in results.toml below.
    ax3 = Axis(fig[1, 2], xticks = (collect(x), labs),
               xlabel = "visc3d_tol [dex of lateral spread]",
               ylabel = "wall clock, full deglaciation [min]")
    kw = .!isnan.(wall)
    barplot!(ax3, collect(x)[kw], wall[kw]; color = C_3D, width = 0.5)

    panels!(fig, [(1, 1), (1, 2)])
    colgap!(fig.layout, 18)

    savefig(fig, "visc3d-tol";
        title = "The lateral-promotion threshold",
        caption = """
        Effect of `visc3d_tol`, the threshold in decades of lateral viscosity spread above
        which a radial element is advanced with laterally-varying rather than radial
        rheology, on the last deglaciation at ℓmax 170 with the 3-D Bagge et al. (2021)
        field. (a) Number of the model's $nt radial elements promoted to the 3-D treatment
        at each threshold; the dashed line is the total element count. (b) Measured wall
        clock for the complete 26 kyr run at each threshold, on 64 threads of one Levante
        node. The production default is 10⁻³.
        Against `tol = 0`, where every Maxwell element is promoted, the approximation
        costs nothing measurable: $(acc_phrase) The threshold is therefore plotted
        against cost alone, because there is no accuracy axis left to plot it against.
        """,
        notes = """
        This is the model-specific approximation that makes the 3-D path affordable, so it
        is the one sensitivity a reviewer is entitled to see quantified rather than
        asserted — and quantified, it turns out to cost nothing. `tol = 0` is the
        reference rather than the default, for the same reason ℓmax 256 is the reference
        in the convergence figure: measuring against the default would make the default
        exact by construction.
        This figure carried an accuracy panel until the degree-1 Z⁶ defect was fixed. It
        showed 0.403 m rms at every threshold, which read as the price of the
        approximation and was nothing of the kind: `tol = 0` promotes 43 elements whose
        viscosity is laterally constant to fourteen digits, so promoting them cannot
        change the answer, and the 0.403 m was the defect those promotions exposed. With it fixed the difference falls to machine precision --
        the two tighter thresholds are bit-identical to `tol = 0` on the model's own
        grid, and the $(@sprintf("%.1e", mx2)) m rms quoted for them here is round-off
        from the regrid the reduction does before differencing -- so the panel was
        removed rather than redrawn.
        Not every radial element can be promoted even at `tol = 0`: elements that are
        elastic rather than Maxwell have no lateral viscosity to vary, which is why (a)
        saturates well below the dashed total.
        10⁻³ and 10⁻² promote the same $(n3[2]) elements and 10⁻¹ promotes $(n3[end]), so
        the threshold is insensitive across two decades. The claim the default has to
        support is the single step from `tol = 0`, and it is
        $(@sprintf("%.2f", wall[1]/wall[2]))× less run time for $(acc_short)
        """,
        source = @__FILE__,
        runs = [Published.blockE_path()])

    for (i, v) in enumerate(e3)
        tag = replace(v.id, "E3-" => "", ".l170" => "")
        record!("blockE", "E3_$(tag)_ne3d" => n3[i])
        isnan(wall[i]) || record!("blockE", "E3_$(tag)_wall_min" => round(wall[i], digits = 1))
        isnan(drsl[i]) || record!("blockE",
            "E3_$(tag)_rms_vs_all_m" => round(drsl[i], digits = 3),
            "E3_$(tag)_uplift_rms_vs_all_mmyr" => round(dupl[i], digits = 4))
    end
    e3
end

# =============================================================================
# Drive, and the table
# =============================================================================
ok1 = fig_lmax()
ok3 = fig_tol()

hhmm(min) = isnan(min) ? "—" : (s = round(Int, 60min); @sprintf("%d:%02d", s ÷ 60, s % 60))

open(joinpath(tables_dir(), "blockE.md"), "w") do io
    println(io, "# Block E — sensitivity\n")
    println(io, "Generated by `analysis/run_blockE.jl`. Both axes are run with the block D ",
                "production configuration otherwise unchanged: GLAC-1D 26 ka → present at ",
                "100 yr, 3-D Bagge (2021) viscosity, `i_eq = 2` against VILMA2's own ",
                "LGM reference state (30 kyr hold), `fe` scheme, rotation on, migrating coastlines.\n")
    println(io, "**Two axes of the original thirteen.** E1 is the numerical convergence of ",
                "the headline result; E3 is `visc3d_tol`, the model's central speed/accuracy ",
                "approximation. E2 (forcing cadence) was dropped because the coupling ",
                "interval is the forcing file's own time axis rather than a namelist ",
                "parameter — see `experiments.md` §E.\n")

    if ok1 !== nothing
        println(io, "## E1 — spectral truncation\n")
        println(io, "Differences are against the most-resolved member, ℓmax ", meta.ref_lmax,
                    ", on the coarser grid of each pair.\n")
        println(io, "| ℓmax | delivered [m] | above-flotation target [m] | delivered/target | rms Δrsl at 21 ka [m] | rms Δuplift [mm/yr] | wall clock |")
        println(io, "|---:|---:|---:|---:|---:|---:|---:|")
        for v in e1
            @printf(io, "| %d | %.2f | %.2f | %.2f %% | %s | %s | %s |\n", v.lmax,
                    v.ocean_mean, v.target, v.closure,
                    isnan(v.rms_rsl) ? "— (ref)" : @sprintf("%.2f", v.rms_rsl),
                    isnan(v.rms_upl) ? "— (ref)" : @sprintf("%.3f", v.rms_upl),
                    hhmm(v.wall_min))
        end
        println(io, "\n**The ℓmax 170 row is about a metre low in both sea-level columns, ",
                    "and that is a reference-state effect rather than a truncation error.** ",
                    "The LGM reference is stored on the ℓmax 170 Gauss grid, so the 170 run ",
                    "reads it natively and every other truncation interpolates it. ",
                    "Interpolation moves the above-flotation volume of the reference ice, ",
                    "so the delivered and target columns move together and the ",
                    "delivered/target ratio stays flat to about a tenth of a percent across ",
                    "the whole ladder. The convergence statement to quote is the ",
                    "64 → 128 → 256 sequence, which shares one interpolation path.\n")
    end

    if ok3 !== nothing
        println(io, "\n## E3 — `visc3d_tol`\n")
        println(io, "Differences are against `tol = 0`, where every Maxwell element is ",
                    "advanced with laterally-varying rheology. The production default is ",
                    "`1e-3`.\n")
        println(io, "| visc3d_tol | elements 3-D | ocean-mean rsl at PD [m] | rms Δrsl at 21 ka [m] | rms Δuplift [mm/yr] | wall clock |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for v in e3
            @printf(io, "| %s | %d | %.2f | %s | %s | %s |\n", v.label,
                    v.ne3d, v.ocean_mean,
                    isnan(v.rms_rsl) ? "— (ref)" : @sprintf("%.3f", v.rms_rsl),
                    isnan(v.rms_upl) ? "— (ref)" : @sprintf("%.4f", v.rms_upl),
                    hhmm(v.wall_min))
        end
    end
end
@info "wrote validation/tables/blockE.md"
