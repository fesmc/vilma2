#!/usr/bin/env julia
#
# run_blockA.jl — parse the Block A benchmark runs staged by
# experiments/run_block_A.sh and produce the validation figures and the summary
# table for the paper's Appendix A / Supporting Information.
#
#   julia --project=. analysis/run_blockA.jl
#
# Figures land in figures/, the summary table in validation/tables/blockA.md.
#
# Reads `data/published/blockA.nc`, not `runs/A`. The benchmark drivers write
# only to stdout, so what is published is the numbers parsed out of it — the
# verdict, the error tables and the SLURM accounting — together with the column
# tables the drivers dump alongside (the profiles and time series behind those
# error numbers, written when VILMA_BENCH_DUMP is set).

using CairoMakie
using Printf
using Statistics

include(joinpath(@__DIR__, "Published.jl"))
include(joinpath(@__DIR__, "FigStyle.jl"))
using .Published
using .FigStyle

const REPO = dirname(@__DIR__)
const SRC  = Published.blockA_path()

runs  = load_blockA()
byid  = Dict(r.id => r for r in runs)
isempty(runs) && error("no published benchmark data; run ./run.sh reduce")

@info "loaded $(length(runs)) run(s) from $(relpath(SRC, REPO))"

# Model against published reference, consistently, throughout every figure.
const C_MODEL = C_1D
const C_REF   = C_3D
# Sequential ramp for the five epochs of the disc transient (light → dark).
const C_TIME  = ["#a8cee8", "#6ba7d4", "#3a7fb8", "#1f5c99", "#0d3b66"]
# The paper's production resolution; drawn on every degree axis for context.
const LMAX_PROD = 170

"Lowest x where y first exceeds tol, or nothing."
function crossing(x, y, tol)
    i = findfirst(>(tol), y)
    i === nothing ? nothing : x[i]
end

# =============================================================================
# A1 — Love numbers over the full degree range (Fig. S1)
# =============================================================================
function fig_love(r)
    el = dumptable(r, "love_elastic.txt")
    fl = dumptable(r, "love_fluid.txt")
    isempty(el) && (@warn "A1: no love_elastic.txt in $(r.id) — rerun with the dump enabled"; return nothing)

    j, h, href, l, lref, k, kref = (el[:, i] for i in 1:7)
    relh = 100 .* abs.(h .- href) ./ abs.(href)
    rell = 100 .* abs.(l .- lref) ./ abs.(lref)
    relk = 100 .* abs.(k .- kref) ./ abs.(kref)

    fig = Figure(size = (980, 720))
    jticks = ([2, 8, 32, 64, 128, 256], ["2", "8", "32", "64", "128", "256"])

    ax1 = Axis(fig[1, 1], xscale = log10, xticks = jticks, xlabel = "harmonic degree j",
               ylabel = "load Love number h", title = "elastic")
    lines!(ax1, j, h, color = C_MODEL, linewidth = 2, label = "VILMA2")
    scatter!(ax1, j[1:8:end], href[1:8:end], color = C_REF, markersize = 6, strokewidth = 0,
             label = "Spada et al. (2011), M3–L70–V01")
    axislegend(ax1, position = :lb)

    ax2 = Axis(fig[1, 2], xscale = log10, xticks = jticks, xlabel = "harmonic degree j",
               ylabel = "load Love number k", title = "elastic")
    lines!(ax2, j, k, color = C_MODEL, linewidth = 2)
    scatter!(ax2, j[1:8:end], kref[1:8:end], color = C_REF, markersize = 6, strokewidth = 0)

    # Error panels on a log axis: the interesting structure spans four decades.
    function errpanel(pos, title, series)
        ax = Axis(fig[pos...], xscale = log10, yscale = log10, xticks = jticks,
                  xlabel = "harmonic degree j", ylabel = "|relative error| [%]", title = title)
        hlines!(ax, [1.0], color = C_GREY, linestyle = :dash)
        vlines!(ax, [LMAX_PROD], color = C_GREY, linestyle = :dot)
        for (name, e, col) in series
            # clamp at the axis floor so exact-zero residuals stay plottable
            lines!(ax, j, max.(e, 1e-4), color = col, linewidth = 2, label = name)
        end
        axislegend(ax, position = :lt)
        ax
    end

    errpanel((2, 1), "elastic",
             [("h", relh, C_MODEL), ("l", rell, "#4c9a2a"), ("k", relk, C_REF)])

    if !isempty(fl)
        jf = fl[:, 1]
        fr = [100 .* abs.(fl[:, 2c] .- fl[:, 2c+1]) ./ abs.(fl[:, 2c+1]) for c in 1:3]
        ax = Axis(fig[2, 2], xscale = log10, yscale = log10, xticks = jticks,
                  xlabel = "harmonic degree j", ylabel = "|relative error| [%]",
                  title = "fluid (relaxed) limit")
        hlines!(ax, [1.0], color = C_GREY, linestyle = :dash)
        vlines!(ax, [LMAX_PROD], color = C_GREY, linestyle = :dot)
        for (name, e, col) in zip(["h", "l", "k"], fr, [C_MODEL, "#4c9a2a", C_REF])
            lines!(ax, jf, max.(e, 1e-4), color = col, linewidth = 2, label = name)
        end
        axislegend(ax, position = :lt)
    end

    xh, xl, xk = crossing(j, relh, 1), crossing(j, rell, 1), crossing(j, relk, 1)
    panels!(fig, [(1, 1), (1, 2), (2, 1), (2, 2)])
    cross(x) = x === nothing ? "never below degree 256" : "at degree $(Int(x))"
    at(lim) = (m = j .<= lim; (maximum(relh[m]), maximum(rell[m]), maximum(relk[m])))
    a170 = at(LMAX_PROD)

    savefig(fig, "love-numbers";
        title = "Load Love numbers against the benchmark table",
        caption = """
        Load Love numbers of the M3–L70–V01 Earth model against the Spada et al. (2011)
        benchmark table, over degrees 2–256. (a, b) The elastic load Love numbers h and k:
        VILMA2 as a line, the benchmark table as points, every eighth degree shown.
        (c) Magnitude of the relative error in the elastic h, l and k. (d) The same in the
        fluid (relaxed) limit. In (c) and (d) the dashed horizontal is 1 % and the dotted
        vertical is ℓmax $LMAX_PROD, the truncation used for the production runs; both
        error axes are logarithmic and residuals are clamped at 10⁻⁴ % so that exact
        zeros remain plottable.
        """,
        notes = """
        This figure is the evidence for a correction the manuscript has to make. The
        model's benchmark documentation claims agreement of about 0.1 % over all degrees
        2–256; that holds for h, which stays at $(@sprintf("%.2f", a170[1])) % to degree
        256, but not for l or k. The elastic error leaves 1 % for l $(cross(xl)) and for k
        $(cross(xk)) — both below the production truncation, where k is off by
        $(@sprintf("%.1f", a170[3])) %. Note that k is small at high degree, so the
        absolute error there is tiny; whether it matters depends on the short-wavelength
        geoid the application needs.
        The test in the model's suite asserts only degrees 2–48, so nothing regressed —
        the claim was simply broader than what was checked. The driver now prints where
        each Love number crosses 1 %.
        """,
        source = @__FILE__,
        runs = [SRC])

    (; ndeg = length(j), cross = (xh, xl, xk), at170 = at(LMAX_PROD), at48 = at(48))
end

# =============================================================================
# A2 — Spada disc load: spatial profiles and the centre transient (Fig. A1)
# =============================================================================
function fig_disc(r)
    els = dumptable(r, "disc_elastic_profile.txt")
    vep = dumptable(r, "disc_ve_profile.txt")
    ser = dumptable(r, "disc_ve_centre_series.txt")
    tab = Published.disc(r).centre
    isempty(vep) && (@warn "A2: no disc_ve_profile.txt in $(r.id)"; return nothing)

    labels = ["0 (elastic)", "1 kyr", "2 kyr", "5 kyr", "10 kyr"]
    θ = vep[:, 1]

    fig = Figure(size = (1000, 760))

    ax1 = Axis(fig[1, 1], xlabel = "colatitude θ [deg]", ylabel = "vertical displacement u [m]",
               title = "relaxation")
    for (i, lab) in enumerate(labels)
        lines!(ax1, θ, vep[:, 1+i], color = C_TIME[i], linewidth = 2, label = lab)
        scatter!(ax1, θ[1:12:end], vep[1:12:end, 1+5+i], color = C_TIME[i], markersize = 6,
                 marker = :circle, strokewidth = 0.5, strokecolor = :white)
    end
    vlines!(ax1, [10.0], color = C_GREY, linestyle = :dash)
    axislegend(ax1, position = :rb)

    # The reference is piecewise linear on a ~1 deg grid (see Published.linear_nodes):
    # only its breakpoints are data. Draw the residual against the full 0.1 deg
    # grid faintly, and mark the residual at the reference's own nodes, which is
    # the part that actually measures the model.
    nodes = [linear_nodes(vep[:, 1+5+i]) for i in 1:5]
    ax2 = Axis(fig[1, 2], xlabel = "colatitude θ [deg]", ylabel = "u − u_ref [m]",
               title = "profile residual")
    hlines!(ax2, [0.0], color = C_GREY, linestyle = :dash)
    for (i, lab) in enumerate(labels)
        resid = vep[:, 1+i] .- vep[:, 1+5+i]
        lines!(ax2, θ, resid, color = (C_TIME[i], 0.30), linewidth = 1.2)
        scatter!(ax2, θ[nodes[i]], resid[nodes[i]], color = C_TIME[i], markersize = 6,
                 strokewidth = 0, label = lab)
    end
    vlines!(ax2, [10.0], color = C_GREY, linestyle = :dash)
    # no legend here: the epochs are the colours of (a), and every corner has data

    if !isempty(els)
        ax3 = Axis(fig[2, 1], xlabel = "colatitude θ [deg]", ylabel = "geoid N [m]",
                   title = "elastic geoid")
        lines!(ax3, els[:, 1], els[:, 4], color = C_MODEL, linewidth = 2, label = "VILMA2")
        scatter!(ax3, els[1:8:end, 1], els[1:8:end, 5], color = C_REF, markersize = 6,
                 strokewidth = 0, label = "Spada et al. (2011)")
        vlines!(ax3, [10.0], color = C_GREY, linestyle = :dash)
        axislegend(ax3, position = :rt)
    end

    ax4 = Axis(fig[2, 2], xlabel = "time [kyr]", ylabel = "centre uplift u(0,t) [m]",
               title = "centre transient")
    if !isempty(ser)
        lines!(ax4, ser[:, 1], ser[:, 2], color = C_MODEL, linewidth = 2,
               label = "VILMA2 (fe, dt = 20 yr)")
    end
    if !isempty(tab)
        scatter!(ax4, tab[:, 1], tab[:, 3], color = C_REF, markersize = 11, strokewidth = 0,
                 label = "Spada et al. (2011)")
    end
    axislegend(ax4, position = :rt)
    panels!(fig, [(1, 1), (1, 2), (2, 1), (2, 2)])

    # Peak-normalized residual per epoch, evaluated at the reference's own data
    # nodes. Comparing on the full 0.1 deg grid instead charges the model for the
    # reference's linear interpolation between nodes; `relall` keeps that figure
    # so the difference between the two measures stays visible.
    #
    # Note also what is and is not asserted upstream: the driver checks the
    # ELASTIC profile spatially and the viscoelastic response only at the centre.
    # The viscoelastic profile away from the centre is newly visible here.
    res    = abs.(vep[:, 2:6] .- vep[:, 7:11])
    peak   = [maximum(abs, vep[:, 1+5+i]) for i in 1:5]
    rel    = [maximum(res[nodes[i], i]) / peak[i] for i in 1:5]
    relall = [maximum(res[:, i]) / peak[i] for i in 1:5]
    worst  = argmax(rel)

    savefig(fig, "disc-benchmark";
        title = "Spada disc load",
        caption = """
        The Spada et al. (2011) disc-load benchmark: a 10° disc of 1000 m of ice applied
        as a Heaviside load on the M3–L70–V01 Earth model, no ocean, at ℓmax 256 with a
        20 yr time step. (a) Vertical displacement against colatitude at five epochs —
        the elastic response and 1, 2, 5 and 10 kyr — as lines for VILMA2 and points
        for the published reference. (b) The residual of those profiles; filled circles
        mark the reference's own data nodes and the faint lines the residual against its
        interpolated grid. (c) The elastic geoid profile. (d) Uplift at the disc centre
        through the relaxation. The dashed vertical in (a)–(c) is the disc edge at 10°.
        """,
        notes = """
        Panel (b) is drawn this way because the published curves are piecewise linear on a
        roughly 1° grid and are distributed interpolated onto a 0.1° grid: only the
        breakpoints carry data. Measured at those $(length(nodes[worst])) nodes the
        peak-normalised residual is
        $(join([@sprintf("%.2f %%", 100r) for r in rel], ", ")) for the five epochs;
        measured against all $(size(vep, 1)) interpolated points it reads
        $(join([@sprintf("%.2f %%", 100r) for r in relall], ", ")). The difference is the
        reference's own interpolation, not model error.
        The model's driver asserts the disc profile only for the elastic epoch and checks
        the viscoelastic response at the centre alone, so the off-centre viscoelastic
        residuals in (b) are a new measurement rather than a re-plot of an asserted one.
        """,
        source = @__FILE__,
        runs = [SRC])

    (; rel, relall, worst, nnode = length(nodes[worst]),
       worstθ = θ[nodes[worst][argmax(res[nodes[worst], worst])]], npoints = size(vep, 1))
end

# =============================================================================
# A3 case A — Martinec (2018) cap loading, no ocean
# =============================================================================
"""
Case A of Martinec et al. (2018): the loading response after 10 kyr against
their fig. 10, from the driver's `martinec_A.txt` dump (colatitude, then
reference and VILMA2 for vertical displacement, horizontal displacement and
geoid).
"""
function fig_martinec_A(r)
    hasdump(r, "martinec_A.txt") ||
        (@warn "A3: no martinec_A.txt in $(r.id) — needs test_benchmark_martinec from 83f63e4 on"; return nothing)
    d = dumptable(r, "martinec_A.txt")
    θ = d[:, 1]
    fields = [("vertical displacement u [m]", 3, 2), ("horizontal displacement [m]", 5, 4),
              ("geoid N [m]", 7, 6)]
    fig = Figure(size = (1250, 620))
    errs = Float64[]
    for (jx, (ylab, cm, cr)) in enumerate(fields)
        ax = Axis(fig[1, jx], ylabel = ylab)
        hidexdecorations!(ax, grid = false, ticks = false)
        lines!(ax, θ, d[:, cm], color = C_MODEL, linewidth = 2, label = "VILMA2")
        stp = max(1, length(θ) ÷ 60)
        scatter!(ax, θ[1:stp:end], d[1:stp:end, cr], color = C_REF, markersize = 6,
                 label = "Martinec et al. (2018)")
        jx == 1 && axislegend(ax, position = :rb)
        peak = maximum(abs.(d[:, cr]))
        res = 100 .* (d[:, cm] .- d[:, cr]) ./ peak
        push!(errs, maximum(abs.(res)))
        axr = Axis(fig[2, jx], xlabel = "colatitude [deg]", ylabel = "residual [% of peak]")
        lines!(axr, θ, res, color = C_MODEL, linewidth = 1.5)
        hlines!(axr, [0.0], color = C_GREY, linestyle = :dash)
        linkxaxes!(ax, axr)
    end
    rowsize!(fig.layout, 2, Relative(0.35))
    panels!(fig, [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3)])
    colgap!(fig.layout, 18); rowgap!(fig.layout, 10)

    savefig(fig, "martinec-A";
        title = "Martinec loading benchmark, case A",
        caption = """
        Case A of the Martinec et al. (2018) benchmark: the response of the M3–L70–V01
        Earth to the spherical-cap ice load after 10 kyr, without ocean, against their
        fig. 10 (the SBK reference). (a) Vertical displacement, (b) horizontal
        displacement and (c) geoid height against colatitude, VILMA2 as the line and the
        reference as points. (d–f) VILMA2 minus the reference, as a percentage of the
        reference peak of each field.
        """,
        notes = """
        The driver asserts these fields against the reference; its maximum errors, as
        printed, are in `tables/blockA.md`. Peak-normalised maxima computed here from the
        dump: $(join([@sprintf("%.2f %%", e) for e in errs], ", ")) for u, horizontal
        displacement and N.
        """,
        source = @__FILE__,
        runs = [SRC])
    (; errs)
end

# =============================================================================
# A3 — Martinec (2018) SLE: the profiles behind the error table (Fig. A2)
# =============================================================================
# fig10/fig12 are constant-longitude cuts (col1 = colatitude); fig11/fig13 are
# constant-colatitude cuts (col1 = 180 + longitude). Same convention as the
# SBK reference files.
const SLE_FIGS = ["fig10", "fig11", "fig12", "fig13"]
sle_xlabel(f) = f in ("fig10", "fig12") ? "colatitude [deg]" : "180 + longitude [deg]"

"""
The per-case profile grid: every Martinec cut, three fields each.

Diagnostic rather than paper material — four of these, one per case, is more
than any appendix needs. The paper figure below carries one case and the error
summary over all of them.
"""
function fig_sle_profiles(id)
    haskey(byid, id) || return nothing
    r = byid[id]
    case = replace(id, "A3." => "")
    have = [f for f in SLE_FIGS if hasdump(r, "sle_$(case)_$(f)_profile.txt")]
    isempty(have) && (@warn "A3: no profile dumps in $id"; return nothing)

    # Columns: 1 col1 | 2,3 u | 4,5 vθ | 6,7 vφ | 8,9 N | 10,11 ss | 12,13 sle
    panels = [("vertical displacement u [m]", 2, 3),
              ("geoid N [m]", 8, 9),
              ("relative sea level [m]", 12, 13)]

    fig = Figure(size = (1100, 260 * length(have)))
    for (i, f) in enumerate(have)
        d = dumptable(r, "sle_$(case)_$(f)_profile.txt")
        for (jx, (ylab, cm, cr)) in enumerate(panels)
            ax = Axis(fig[i, jx], xlabel = sle_xlabel(f), ylabel = ylab,
                      title = jx == 1 ? "$f" : "")
            lines!(ax, d[:, 1], d[:, cm], color = C_MODEL, linewidth = 2, label = "VILMA2")
            scatter!(ax, d[1:8:end, 1], d[1:8:end, cr], color = C_REF, markersize = 6,
                     strokewidth = 0, label = "Martinec et al. (2018)")
            (i == 1 && jx == 1) && axislegend(ax, position = :rb)
        end
    end
    panels!(fig, [(i, jx) for i in 1:length(have) for jx in 1:3])

    savefig(fig, "sle-profiles-$case";
        title = "Martinec SLE benchmark, case $case — all cuts",
        caption = """
        Case $case of the Martinec et al. (2018) sea-level-equation benchmark, on every cut
        the reference provides ($(join(have, ", "))), against the SBK reference curves.
        Columns are vertical displacement, geoid height and relative sea level;
        VILMA2 is the line and the reference the points. Cuts fig10 and fig12 are at
        constant longitude and are plotted against colatitude; fig11 and fig13 are at
        constant colatitude and are plotted against 180 + longitude, following the
        reference files' own convention.
        """,
        notes = """
        A diagnostic: the paper's benchmark figure carries one case and the error summary
        over all four, which is what an appendix needs. This grid exists so that a
        disagreement can be located on a specific cut and field rather than inferred from
        a bar chart.
        """,
        source = @__FILE__,
        runs = [SRC])
end

"""
The A3 paper figure: one case in profile, the run history for all cases, and the
error summary that the other three cases are represented by.
"""
function fig_sle_benchmark(cases; show_case = "D3")
    id = "A3.$show_case"
    haskey(byid, id) || (@warn "A3: $id not available for the benchmark figure"; return nothing)
    r = byid[id]
    cut = first([f for f in SLE_FIGS if hasdump(r, "sle_$(show_case)_$(f)_profile.txt")], 1)
    isempty(cut) && (@warn "A3: no profile dump for $id"; return nothing)
    d = dumptable(r, "sle_$(show_case)_$(cut[1])_profile.txt")

    hist = [(replace(c, "A3." => ""), dumptable(byid[c], "sle_$(replace(c, "A3." => ""))_history.txt"))
            for c in cases if haskey(byid, c)]
    filter!(x -> !isempty(x[2]), hist)
    errs = [(replace(c, "A3." => ""), Published.sle(byid[c])) for c in cases if haskey(byid, c)]
    filter!(x -> !isempty(x[2].table), errs)

    fig = Figure(size = (1250, 720))
    for (jx, (ylab, cm, cr)) in enumerate([("vertical displacement u [m]", 2, 3),
                                           ("geoid N [m]", 8, 9),
                                           ("relative sea level [m]", 12, 13)])
        ax = Axis(fig[1, jx], xlabel = sle_xlabel(cut[1]), ylabel = ylab)
        lines!(ax, d[:, 1], d[:, cm], color = C_MODEL, linewidth = 2, label = "VILMA2")
        scatter!(ax, d[1:8:end, 1], d[1:8:end, cr], color = C_REF, markersize = 6,
                 strokewidth = 0, label = "Martinec et al. (2018)")
        jx == 1 && axislegend(ax, position = :rb)
    end

    ax4 = Axis(fig[2, 1], xlabel = "time [kyr]", ylabel = "eustatic sea level [m]")
    # The per-step residual sits at double-precision noise (1e-16 and exact
    # zeros), so plotting it raw is a hairball that says nothing. The running
    # maximum answers the question that matters -- does it ever drift? -- as a
    # clean monotone staircase.
    ax5 = Axis(fig[2, 2], xlabel = "time [kyr]", yscale = log10,
               ylabel = "running max |ocean-mass residual|")
    # E2 and F1 end within 0.4 m of each other, so each case also gets its own
    # dash pattern and the later ones are drawn thinner, on top.
    styles = [:solid, :dash, :dot, :dashdot]
    for (i, (case, h)) in enumerate(hist)
        col = Makie.wong_colors()[i]
        lines!(ax4, h[:, 2], h[:, 4], color = col, linewidth = 3.0 - 0.5 * (i - 1),
               linestyle = styles[mod1(i, length(styles))], label = case)
        # Floor at 1e-18: the first steps are exactly zero, which a log axis
        # would otherwise stretch to meaninglessness.
        lines!(ax5, h[:, 2], max.(accumulate(max, abs.(h[:, 5])), 1e-18), color = col,
               linewidth = 2, label = case)
    end
    ylims!(ax5, 1e-18, 1e-13)
    hlines!(ax5, [eps(Float64)], color = C_GREY, linestyle = :dash)
    axislegend(ax4, position = :rc)

    ax6 = Axis(fig[2, 3], ylabel = "peak-normalized max error [%]",
               xticks = (1:length(errs), [c for (c, _) in errs]))
    fields = ["u", "N", "sle"]
    fcols  = [2, 5, 7]                      # columns of Published.sle's table
    for (c, (name, fc)) in enumerate(zip(fields, fcols))
        xs = (1:length(errs)) .+ (c - (length(fields) + 1) / 2) * (0.8 / length(fields))
        ys = [maximum(dd.table[:, fc]) for (_, dd) in errs]
        barplot!(ax6, xs, ys, width = 0.8 / length(fields),
                 color = Makie.wong_colors()[c], label = name)
    end
    axislegend(ax6, position = :lt)

    panels!(fig, [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3)])
    colgap!(fig.layout, 18); rowgap!(fig.layout, 14)

    savefig(fig, "sle-benchmark";
        title = "Martinec sea-level-equation benchmark",
        caption = """
        The Martinec et al. (2018) sea-level-equation benchmark, at ℓmax 128 with a 20 yr
        step and the subgrid water-column ocean load. (a–c) Case $show_case on the
        $(cut[1]) cut: vertical displacement, geoid height and relative sea level, with
        VILMA2 as the line and the SBK reference as points. (d) Eustatic sea level
        through the run for all four cases. (e) Running maximum of the absolute
        ocean-mass residual for the same runs; the dashed horizontal is double-precision
        machine epsilon. (f) Peak-normalised maximum error against the reference curves,
        per case, for displacement, geoid and relative sea level.
        """,
        notes = """
        Panel (e) is a running maximum rather than the per-step residual because the
        per-step value sits at double-precision noise — exact zeros and values near
        10⁻¹⁶ — where a raw trace carries no information. What the panel is asked to show
        is whether the residual ever drifts, and a monotone staircase answers that
        directly.
        Case $show_case is shown in profile because one case is what an appendix figure
        can carry legibly; (f) is how the other three are represented, and the full
        per-case profile grids are in `figures/diagnostics/sle-profiles-*`.
        The reference curves are tabulated rather than interpolated, so unlike the disc
        benchmark these residuals are not limited by the reference's own resolution.
        """,
        source = @__FILE__,
        runs = [SRC])

    Dict("A3." * c => (; maxerr_u = maximum(dd.table[:, 2]), maxerr_sle = maximum(dd.table[:, 7]),
                         esl = isempty(dd.eustatic) ? NaN : dd.eustatic[1]) for (c, dd) in errs)
end

# =============================================================================
# A4 — rotation / polar motion (Fig. S4)
# =============================================================================
function fig_rotation(r)
    d = Published.rotation(r)
    bodies = sort(collect(keys(d.bodies)))
    series = Dict(b => dumptable(r, "rotation_$(b)_series.txt") for b in bodies)
    isempty(bodies) && (@warn "A4: no polar-motion table in $(r.id)"; return nothing)

    fig = Figure(size = (1050, 380))
    ax1 = Axis(fig[1, 1], xlabel = "time [kyr]", ylabel = "|m(t)| [deg]",
               title = "amplitude")
    # m1 and m2 grow almost proportionally here, so an m1-m2 path is a
    # near-straight line with a degenerate aspect ratio; show the two
    # components against time instead, which is what carries the information.
    ax2 = Axis(fig[1, 2], xlabel = "time [kyr]", ylabel = "polar motion component [deg]",
               title = "components")
    ax3 = Axis(fig[1, 3], xlabel = "time [kyr]", ylabel = "relative error [%]",
               title = "error against the reference")
    hlines!(ax3, [0.0], color = C_GREY, linestyle = :dash)

    for (i, b) in enumerate(bodies)
        col = i == 1 ? C_MODEL : "#4c9a2a"
        s = series[b]
        if !isempty(s)
            lines!(ax1, s[:, 1], s[:, 4], color = col, linewidth = 2, label = "$b (model)")
            lines!(ax2, s[:, 1], s[:, 2], color = col, linewidth = 2, label = "$b  m₁")
            lines!(ax2, s[:, 1], s[:, 3], color = col, linewidth = 2, linestyle = :dash,
                   label = "$b  m₂")
        end
        tab = d.bodies[b]
        isempty(tab) && continue
        scatter!(ax1, tab[:, 1], tab[:, 3], color = C_REF, markersize = 9,
                 marker = i == 1 ? :circle : :rect, strokewidth = 0,
                 label = "$b (Spada Table 14)")
        # Recomputed from the tabulated values: the driver's printed rel.err
        # column carries only three decimals, i.e. 0.1 % resolution.
        rel = 100 .* (tab[:, 2] .- tab[:, 3]) ./ abs.(tab[:, 3])
        scatterlines!(ax3, tab[:, 1], rel, color = col, markersize = 8, label = b)
    end
    axislegend(ax1, position = :rb)
    axislegend(ax2, position = :rc)      # between the m1 and m2 curves
    axislegend(ax3, position = :rt)
    panels!(fig, [(1, 1), (1, 2), (1, 3)])

    maxrel = maximum(maximum(abs, 100 .* (t[:, 2] .- t[:, 3]) ./ abs.(t[:, 3]))
                     for t in values(d.bodies) if !isempty(t))

    savefig(fig, "polar-motion";
        title = "Rotational feedback",
        caption = """
        Polar motion driven by surface loading, against test 3/2 of Spada et al. (2011),
        with rotation on and the sea-level equation off. (a) Amplitude |m(t)| for each
        load geometry, model as lines and the published Table 14 values as points.
        (b) The two components m₁ (solid) and m₂ (dashed) against time. (c) Relative
        error against the tabulated values, recomputed from the table rather than read
        from the driver's printed column, which carries only three decimals. The maximum
        is $(@sprintf("%.2f", maxrel)) %. The degree-2 tidal Love numbers this run uses
        are elastic kᵀₑ = $(@sprintf("%.4f", d.kTe)) and secular k_s =
        $(@sprintf("%.4f", d.k_s)).
        """,
        notes = """
        Panel (b) shows the two components against time rather than the m₁–m₂ path
        because in this configuration m₁ and m₂ grow almost proportionally: the path is a
        near-straight line with a degenerate aspect ratio and carries less than the two
        traces do.
        This is the rotation test without the sea-level equation. The coupled case — the
        rotational feedback entering the SLE — is a separate run in the suite (A4b) and
        is reported in the summary table as a pass/fail with its mass residual rather
        than as a figure.
        """,
        source = @__FILE__,
        runs = [SRC])

    (; bodies, maxrel)
end

# =============================================================================
# A5 — low-viscosity zone (Fig. S5)
# =============================================================================
function fig_lvz(r)
    d = Published.lvz(r)
    cases = [("homog", "homogeneous η = 10²¹ Pa s", C_MODEL),
             ("lvz",   "with LVZ (10¹⁹, 70–170 km)", C_REF)]
    prof = Dict(c => dumptable(r, "lvz_$(c)_profile.txt") for (c, _, _) in cases)
    ser  = Dict(c => dumptable(r, "lvz_$(c)_series.txt")  for (c, _, _) in cases)
    all(isempty, values(prof)) && (@warn "A5: no LVZ dumps in $(r.id)"; return nothing)

    fig = Figure(size = (1000, 400))

    ax1 = Axis(fig[1, 1], xlabel = "distance from load centre x [km]",
               ylabel = "vertical displacement u [m]",
               title = "profile at 200 yr")
    for (c, lab, col) in cases
        p = prof[c]
        isempty(p) && continue
        lines!(ax1, p[:, 1], p[:, 2], color = col, linewidth = 2, label = lab)
    end
    vlines!(ax1, [100.0], color = C_GREY, linestyle = :dash)
    axislegend(ax1, position = :rb)

    ax2 = Axis(fig[1, 2], xlabel = "time [yr]", ylabel = "centre uplift u(0,t) [m]",
               title = "centre transient")
    for (c, lab, col) in cases
        s = ser[c]
        isempty(s) && continue
        lines!(ax2, s[:, 1], s[:, 2], color = col, linewidth = 2, label = lab)
    end
    for (c, _, col) in cases
        key = c == "homog" ? "homogeneous" : "with LVZ"
        haskey(d.cases, key) || continue
        scatter!(ax2, [200.0], [d.cases[key][2]], color = col, markersize = 11,
                 marker = :star5, strokewidth = 0)
    end
    vlines!(ax2, [100.0], color = C_GREY, linestyle = :dash)
    axislegend(ax2, position = :lb)
    panels!(fig, [(1, 1), (1, 2)])

    amp = isempty(d.amplification) ? "not reported" :
          @sprintf("%.3f against a published %.2f", d.amplification[1], d.amplification[2])
    errs = join([@sprintf("%s %.4f m against %.2f m (%.1f %%)", n, v[1], v[2], v[3])
                 for (n, v) in sort(collect(d.cases), by = first)], "; ")

    savefig(fig, "lvz-benchmark";
        title = "Low-viscosity zone",
        caption = """
        The low-viscosity-zone benchmark of Weerdesteijn et al. (2023, §5.2) at ℓmax 512:
        a 100 km disc loaded from 0 to 100 m over 100 yr and then held to 200 yr, solved
        on a homogeneous Earth (η = 10²¹ Pa s) and on one carrying a low-viscosity zone
        (10¹⁹ Pa s between 70 and 170 km depth). (a) Vertical displacement against
        distance from the load centre at 200 yr. (b) Uplift at the load centre against
        time; stars are the published reference values at 200 yr. The dashed vertical in
        (a) is the disc edge and in (b) the moment the load stops growing. Centre
        displacements are $errs, and the amplification U_LVZ/U_homog is $amp.
        """,
        notes = """
        This is the only published-benchmark validation of the laterally-varying-viscosity
        path, and therefore of the capability the paper is built on. It had a build target
        but was in no test target until model `1911b4f`; it is now part of
        `make check-slow`, which is why the block A stage submits that suite to the
        compute queue rather than running it on a login node.
        """,
        source = @__FILE__,
        runs = [SRC])

    (; cases = d.cases, amplification = d.amplification)
end

# =============================================================================
# Summary table
# =============================================================================
function write_summary(runs, extras)
    hhmmss(s) = s === nothing ? "—" : @sprintf("%02d:%02d:%02d", s ÷ 3600, (s % 3600) ÷ 60, s % 60)
    io = IOBuffer()

    println(io, "# Block A — community benchmarks, run status\n")
    println(io, "Generated by `analysis/run_blockA.jl` from `", relpath(SRC, REPO), "`.")
    println(io)
    println(io, "| Run | Verdict | Wall | Peak RSS | Node | Job |")
    println(io, "|---|---|---:|---:|---|---:|")
    for r in runs
        m = r.slurm
        rss = haskey(m, "maxrss_kb") ? @sprintf("%.1f MB", m["maxrss_kb"] / 1024) : "—"
        println(io, "| `", r.id, "` | ", verdict(r), " | ",
                hhmmss(get(m, "elapsed_s", nothing)), " | ", rss, " | ",
                get(m, "node", "—"), " | ", get(m, "jobid", "—"), " |")
    end

    println(io, "\n**", count(passed, runs), " of ", length(runs), " runs passed.**\n")

    if !isempty(extras)
        println(io, "## Key quantitative results\n")
        for (k, v) in extras
            println(io, "- **", k, "** — ", v)
        end
    end

    path = joinpath(tables_dir(), "blockA.md")
    write(path, String(take!(io)))
    @info "wrote $path"
end

# =============================================================================
# Driver
# =============================================================================
extras = Pair{String,String}[]
fmt(x) = x === nothing ? "never (≤256)" : "j=$(Int(x))"

if haskey(byid, "A1.love")
    s = fig_love(byid["A1.love"])
    if s !== nothing
        push!(extras, "A1 Love numbers" =>
            @sprintf("degrees 2–%d. Max |rel. error| for j≤48: h %.3f %%, l %.3f %%, k %.3f %%. At the production ℓmax=%d: h %.3f %%, l %.2f %%, k %.2f %%. Leaves 1 %% at h %s, l %s, k %s.",
                     s.ndeg + 1, s.at48..., LMAX_PROD, s.at170..., fmt.(s.cross)...))
    end
end

if haskey(byid, "A2.disc")
    s = fig_disc(byid["A2.disc"])
    if s !== nothing
        labs = ["elastic", "1 kyr", "2 kyr", "5 kyr", "10 kyr"]
        push!(extras, "A2 disc load" =>
            @sprintf("peak-normalized max residual per epoch, **at the reference's own ~%d data nodes**: %s (worst %s at θ=%.1f°). Evaluated instead against all %d interpolated points it reads %s — the difference is the reference's linear interpolation, not model error. The driver asserts the profile only for the elastic epoch; the viscoelastic profile is checked upstream at the centre alone.",
                     s.nnode,
                     join([@sprintf("%s %.2f %%", l, 100r) for (l, r) in zip(labs, s.rel)], ", "),
                     labs[s.worst], s.worstθ, s.npoints,
                     join([@sprintf("%.2f %%", 100r) for r in s.relall], "/")))
    end
end

if haskey(byid, "A3.caseA")
    fig_martinec_A(byid["A3.caseA"])
    d = Published.martinec_caseA(byid["A3.caseA"])
    !isempty(d.fields) && push!(extras, "A3 case A" =>
        join([@sprintf("%s %.2f %%", f, e) for (f, e) in zip(d.fields, d.errs)], ", "))
end

let cases = ["A3.C2", "A3.D3", "A3.E2", "A3.F1"]
    s = fig_sle_benchmark(cases)
    for c in cases
        fig_sle_profiles(c)
    end
    if s !== nothing
        for (id, v) in sort(collect(s), by = first)
            push!(extras, "A3 $(replace(id, "A3." => ""))" =>
                @sprintf("max error u %.2f %%, sle %.2f %%; esl %.2f m", v.maxerr_u, v.maxerr_sle, v.esl))
        end
    end
end

if haskey(byid, "A4.rotation")
    s = fig_rotation(byid["A4.rotation"])
    s !== nothing && push!(extras, "A4 polar motion" =>
        @sprintf("%s; max |rel. error| in |m(t)| %.2f %%", join(s.bodies, " + "), s.maxrel))
end

if haskey(byid, "A5.lvz")
    s = fig_lvz(byid["A5.lvz"])
    s !== nothing && push!(extras, "A5 LVZ" =>
        join([@sprintf("%s %.4f m (ref %.2f m, %.1f %%)", n, v[1], v[2], v[3])
              for (n, v) in sort(collect(s.cases), by = first)], "; "))
end

write_summary(runs, extras)
@info "Block A analysis complete — $(count(passed, runs))/$(length(runs)) passed"
