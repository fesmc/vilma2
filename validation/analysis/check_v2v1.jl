#!/usr/bin/env julia
#
# check_v2v1.jl — VILMA2 against VILMA1 on the last deglaciation
#.
#
#   julia --project=. analysis/check_v2v1.jl
#
# TWO PAIRS -- the pruned list. Both are referenced to VILMA2's own
# LGM reference state (i_eq = 2; the 26 ka load held 30 kyr) and both run with the degree-1 frames
# matched (deg1_frame = "cm").
#
# Why i_eq = 2 and nothing else: with an LGM datum VILMA1 begins at its own
# reference with no first-step jump, and VILMA2 begins at an equilibrium of
# its own Earth model. Neither needs a spin-up, so the initialization asymmetry
# that dogged the i_eq = 1 comparison -- V2 entering spun up via I1 while V1
# absorbed the whole LGM-minus-present ice at t0 -- does not arise. That datum
# only became usable once the grounded-ice mask defect was fixed (model 92424dc).
#
# Why "cm": the default "cf" leaves rsl with no degree 1 at all, while VILMA1
# reports geocenter motion. On the disc benchmark that convention difference was
# 98-99.7 % of the V2-V1 residual.
#
# The seven other pairs this table used to carry -- pd, lgm, uniform, matched,
# matched2 and their .cm variants -- were the datum-hunting exercise that the
# mask fix retired. They are gone rather than left skipping, so nobody re-runs
# them looking for a discrepancy that no longer exists.
#
# WHAT IS AND IS NOT LIKE-FOR-LIKE
#
#   same:  ice forcing (GLAC-1D, remapped by the same driver), reference state
#          (one file, read by both), output grid (model Gauss, lmax 170),
#          degree-1 frame, ocean-function and bsl diagnostics (computed by
#          VILMA2 for BOTH backends, by the same formula), and -- since
#          experiments/make_vilma1_visc.jl -- the viscosity field, clamped
#          identically. VILMA1 reads its NetCDF directly, so that last one had to
#          be built, not assumed.
#
#   NOT the same:
#     * LOAD PATH WITHIN A COUPLING INTERVAL. VILMA2 ramps the load
#       linearly across the 100 yr interval and sub-steps adaptively; VILMA1
#       holds it at the interval value and takes vilma_v1_nsub uniform sub-steps
#       (12 here). This is a documented asymmetry, not a bug in either code.
#     * RADIAL DISCRETISATION. Both integrate the same viscosity field, but onto
#       their own element grids.

using CairoMakie
using Printf
using Statistics

include(joinpath(@__DIR__, "Published.jl"))
include(joinpath(@__DIR__, "FigStyle.jl"))
using .Published
using .FigStyle

const REPO = dirname(@__DIR__)

const PAIRS = [
    (tag = "3d",  slug = "v2v1-maps-3d",  v2 = RUN_V2_3D,  v1 = RUN_V1_3D,
     lab = "3-D Bagge (2021) viscosity"),
    (tag = "1db", slug = "v2v1-maps-1db", v2 = RUN_V2_1DB, v1 = RUN_V1_1DB,
     lab = "1Db, the radial mean of that field"),
    # V7's control: the 3-D pair with VILMA2's toroidal coupling off. A
    # diagnostic (its slug is not in DOCS_FIGURES), and a third line in the summary.
    (tag = "3d_notor", slug = "v2v1-maps-3d-notor", v2 = RUN_V2_3D_NOTOR, v1 = RUN_V1_3D,
     lab = "3-D Bagge (2021) viscosity, without VILMA2's toroidal coupling"),
]

# How each pair is drawn in the summary figure: colour and legend label.
const PAIR_STYLE = Dict("3d"       => (col = C_3D,  lab = "3-D"),
                        "1db"      => (col = C_1D,  lab = "1Db"),
                        "3d_notor" => (col = C_3D_DARK, lab = "3-D, no toroidal"))

# V7: the toroidal effect, tor − notor, against the residual it is meant to close.
# The production 3-D run carries the toroidal field; the control is the same
# head, reference and namelist but for l_toroidal.
const TOR, NOTOR, V3D = RUN_V2_3D, RUN_V2_3D_NOTOR, RUN_V1_3D

# Every mark must be one of `Published.EPOCHS`: the table is evaluated from the
# published 2-D fields, and Published.slice_at refuses an epoch that was not
# reduced rather than silently returning its neighbour.
const T_MARKS = [-21000.0, -14000.0, -8000.0, -2000.0]

stats(d) = (; rms = sqrt(mean(d .^ 2)), mx = maximum(abs.(d)), p99 = quantile(abs.(vec(d)), 0.99))

# `deg01` — the area-weighted degree-0/1 fit — is in Published rather than here.
# It is what decides whether the V2−V1 residual is a solver difference or a frame
# convention difference, so it travels with the published data for a reader
# re-deriving the residual offline. See its docstring in analysis/Published.jl.

# =============================================================================
# The map figure, one per pair
# =============================================================================
function maps_figure(p, f, v, SC)
    dupl = uplift_rate(f) .- uplift_rate(v)
    # All four limits come from FigStyle.deglac_scales, computed over the whole
    # block D matrix. Previously each pair autoscaled itself, so the 3-D and 1Db
    # versions of this figure were on DIFFERENT scales (±151 m against ±115 m)
    # and neither could be weighed against the lateral-signal figure. The signal
    # and difference ranges are still separate -- on the V2/V1 range a residual of
    # a few per cent renders as an empty white panel, and the figure would
    # assert agreement it had not shown -- but each is now one range for the set.
    rlim,  ulim  = SC.rsl,  SC.upl
    dlim,  dulim = SC.drsl, SC.dupl

    fig = Figure(size = (1450, 830))
    mkax(pos, title; showx = false, showy = false) =
        (a = Axis(fig[pos...]; title, titlealign = :center, aspect = DataAspect(),
                  xticks = ([-90, 0, 90], ["90°W", "0°", "90°E"]),
                  yticks = ([-60, 0, 60], ["60°S", "0°", "60°N"]),
                  xticklabelsize = TICK_FS, yticklabelsize = TICK_FS,
                  xticklabelsvisible = showx, yticklabelsvisible = showy);
         limits!(a, -180, 180, -90, 90); a)

    hm_rsl = hm_upl = hm_drsl = hm_dupl = nothing
    for (col, (r, name)) in enumerate(((f, "VILMA2"), (v, "VILMA1")))
        for (row, t) in enumerate((-21000.0, -8000.0))
            a = mkax((row, col), row == 1 ? name : ""; showy = col == 1)
            hm_rsl = heatmap!(a, r.lon, r.lat, rsl_anomaly(r, t);
                              colormap = CMAP_DIV, colorrange = rlim)
        end
        a = mkax((3, col), ""; showx = true, showy = col == 1)
        hm_upl = heatmap!(a, r.lon, r.lat, uplift_rate(r); colormap = :curl, colorrange = ulim)
    end
    for (row, t) in enumerate((-21000.0, -8000.0))
        a = mkax((row, 4), row == 1 ? "VILMA2 − VILMA1" : "")
        hm_drsl = heatmap!(a, f.lon, f.lat, rsl_anomaly(f, t) .- rsl_anomaly(v, t);
                           colormap = CMAP_DIV, colorrange = dlim)
    end
    a = mkax((3, 4), ""; showx = true)
    hm_dupl = heatmap!(a, f.lon, f.lat, dupl; colormap = :curl, colorrange = dulim)

    Colorbar(fig[1:2, 3], hm_rsl,  label = "rsl vs present day [m]")
    Colorbar(fig[3, 3],   hm_upl,  label = "uplift rate [mm/yr]")
    Colorbar(fig[1:2, 5], hm_drsl, label = "Δ rsl [m]")
    Colorbar(fig[3, 5],   hm_dupl, label = "Δ uplift rate [mm/yr]")

    # Row labels: the epoch, once, on the left of the row.
    for (row, lab) in enumerate(("21 ka", "8 ka", "present day,\nuplift rate"))
        Label(fig[row, 0], lab, rotation = π/2, fontsize = 13, tellheight = false)
    end
    panels!(fig, [(r, c) for r in 1:3 for c in (1, 2, 4)])
    rowgap!(fig.layout, 6); colgap!(fig.layout, 6)

    savefig(fig, p.slug;
        title = "VILMA2 against VILMA1 — $(p.lab)",
        caption = """
        The last deglaciation solved twice on the identical problem — same GLAC-1D
        forcing, same reference state file, same $(p.lab), same ℓmax 170 grid, same
        degree-1 frame — by VILMA2 (left column) and by VILMA1 (middle column),
        with their difference at right. (a, d, g) VILMA2; (b, e, h) VILMA1;
        (c, f, i) VILMA2 minus VILMA1. Top and middle rows are relative sea level
        at 21 ka and 8 ka relative to each run's own present day; the bottom row is the
        present-day bedrock uplift rate. The two model columns share a colour scale;
        the difference column has its own, narrower one, given by its own colour bar.
        """,
        notes = """
        The difference column is on a separate scale on purpose: at the V2/V1 scale the
        residual renders as a blank white panel, which would show agreement rather than
        demonstrate it. Read the scale bars before comparing panels across columns.
        Two things are deliberately not matched and account for part of the residual:
        VILMA2 ramps the ice load linearly across each 100 yr coupling interval and
        sub-steps adaptively while VILMA1 holds the load and takes 12 uniform sub-steps,
        and the two codes discretise the same viscosity field onto their own radial
        element grids. The `C_ocean` and `bsl` diagnostics are computed by VILMA2
        for both backends, so they are not independent between the columns.
        """,
        source = @__FILE__,
        runs = [f.path, v.path])
end

# =============================================================================
# The summary figure: how big the disagreement is, and what it is made of
# =============================================================================
function summary_figure(series)
    # The degree-0/1 bound the notes quote, measured over the series drawn.
    d01 = maximum(max(maximum(abs.(s.deg0)), maximum(s.deg1)) for (_, s) in series)

    fig = Figure(size = (1150, 720))

    ax1 = Axis(fig[1, 1], xlabel = "time [ka]", ylabel = "rms |V2 − V1| [m]")
    ax2 = Axis(fig[1, 2], xlabel = "time [ka]", ylabel = "rms |V2 − V1| / max |V2| [%]")
    ax3 = Axis(fig[2, 1], xlabel = "time [ka]", ylabel = "|V2 − V1| content [m]")
    ax4 = Axis(fig[2, 2], xlabel = "time [ka]", ylabel = "ocean-function mismatch [% of grid]")

    for (tag, s) in series
        col, lab = PAIR_STYLE[tag]
        t = s.time ./ 1000
        lines!(ax1, t, s.rms;        color = col, label = lab)
        lines!(ax2, t, 100 .* s.rel; color = col, label = lab)
        lines!(ax3, t, abs.(s.deg0); color = col, linestyle = :solid,   label = "$lab, |degree 0|")
        lines!(ax3, t, s.deg1;       color = col, linestyle = :dashdot, label = "$lab, |degree 1|")
        lines!(ax4, t, s.cmis;       color = col, label = lab)
    end
    axislegend(ax1, position = :rt)
    axislegend(ax3, position = :rt)
    linkxaxes!(ax1, ax2, ax3, ax4)
    panels!(fig, [(1, 1), (1, 2), (2, 1), (2, 2)])
    colgap!(fig.layout, 18); rowgap!(fig.layout, 14)

    savefig(fig, "v2v1-summary";
        title = "Size and composition of the VILMA2 − VILMA1 residual",
        caption = """
        The VILMA2 minus VILMA1 residual in relative sea level through the last
        deglaciation, for $(join([lowercase(PAIR_STYLE[tag].lab) == "1db" ? "1Db, the radial mean of the 3-D field (blue)" :
                                  tag == "3d" ? "the 3-D viscosity field (red)" :
                                  "the 3-D field with VILMA2's toroidal coupling off (dark red)"
                                  for (tag, _) in series], ", ")), all at ℓmax 170 on the
        identical reference state. (a) Root-mean-square
        residual over the global grid. (b) The same, as a percentage of the maximum
        |relative sea level| that VILMA2 produces at that epoch. (c) The degree-0
        (solid) and degree-1 (dash-dotted) content of the residual, from an area-weighted
        least-squares fit. (d) Fraction of grid cells on which the two runs' ocean
        functions disagree.
        """,
        notes = """
        Panel (c) is the check that the residual is a solver difference rather than a
        convention difference. With VILMA2 in its `cf` frame and VILMA1 reporting
        geocentre motion, degree 1 was about half the deglacial residual; with both
        in the `cm` frame degrees 0 and 1 stay below
        $(@sprintf("%.2f", d01)) m in every pair, against a peak rms of
        $(join([@sprintf("%.1f m (%s)", maximum(s.rms), PAIR_STYLE[tag].lab)
                for (tag, s) in series], ", ")). The degree-2, order-1 pattern of the
        rotational feedback that dominated earlier passes is gone since the rotation
        fixes of 304f4cb (`tables/v2v1.md`).
        Sea level is referenced to each run's own present day, so the residual at 0 ka
        is zero by construction and (b) is undefined there.
        """,
        source = @__FILE__,
        runs = [Published.v2v1_path()])
end

# =============================================================================
# V7: the toroidal effect, against the residual it is meant to close
# =============================================================================
"""
    toroidal_figure(tr, nt, v, SC)

Three columns: VILMA2 without the toroidal coupling minus VILMA1, the same
with it, and the toroidal effect itself (with − without). Rows are relative sea
level at 21, 14 and 8 ka and the present-day uplift rate. Every panel is on the
shared block D difference range, so the three columns weigh against each other
and against the lateral-signal figure.
"""
function toroidal_figure(tr, nt, v, SC)
    cols = [("VILMA2 − VILMA1,\nno toroidal coupling", nt, v),
            ("VILMA2 − VILMA1,\ntoroidal coupling", tr, v),
            ("toroidal coupling:\nwith − without", tr, nt)]
    rows = [("21 ka", r -> rsl_anomaly(r, -21000.0)),
            ("14 ka", r -> rsl_anomaly(r, -14000.0)),
            ("8 ka",  r -> rsl_anomaly(r, -8000.0)),
            ("present day,\nuplift rate", r -> uplift_rate(r))]
    rmsof(x) = sqrt(mean(x .^ 2))

    fig = Figure(size = (1400, 1050))
    hm_r = hm_u = nothing
    for (i, (rl, get)) in enumerate(rows), (j, (cl, a, b)) in enumerate(cols)
        d = get(a) .- get(b)
        isupl = i == length(rows)
        ax = Axis(fig[i, j]; title = @sprintf("%srms %.3g %s", i == 1 ? cl * "\n" : "",
                                              rmsof(d), isupl ? "mm/yr" : "m"),
                  titlesize = 12, titlealign = :center, aspect = DataAspect(),
                  xticks = ([-90, 0, 90], ["90°W", "0°", "90°E"]),
                  yticks = ([-60, 0, 60], ["60°S", "0°", "60°N"]),
                  xticklabelsize = TICK_FS, yticklabelsize = TICK_FS,
                  xticklabelsvisible = isupl, yticklabelsvisible = j == 1)
        limits!(ax, -180, 180, -90, 90)
        h = heatmap!(ax, a.lon, a.lat, d; colormap = isupl ? :curl : CMAP_DIV,
                     colorrange = isupl ? SC.dupl : SC.drsl)
        isupl ? (hm_u = h) : (hm_r = h)
    end
    for (i, (rl, _)) in enumerate(rows)
        Label(fig[i, 0], rl, rotation = π/2, fontsize = 13, tellheight = false)
    end
    Colorbar(fig[1:3, 4], hm_r, label = "Δ rsl [m]")
    Colorbar(fig[4, 4],   hm_u, label = "Δ uplift rate [mm/yr]")
    panels!(fig, [(r, c) for r in 1:length(rows) for c in 1:3])
    rowgap!(fig.layout, 6); colgap!(fig.layout, 6)

    savefig(fig, "toroidal-effect";
        title = "The toroidal coupling against the VILMA2 − VILMA1 residual",
        caption = """
        The 3-D Bagge (2021) deglaciation at ℓmax 170, run by VILMA2 with and
        without its toroidal coupling from one binary, against the matched VILMA1 run.
        Left: VILMA2 without the toroidal coupling minus VILMA1. Middle: the same
        with it. Right: the toroidal effect, with minus without. Top three rows are
        relative sea level relative to each run's own present day at 21, 14 and 8 ka;
        the bottom row is the present-day bedrock uplift rate. All panels share the
        block D difference range.
        """,
        notes = """
        The run with the toroidal coupling is the production 3-D run, so the middle
        column is the 3-D pair of the V2−V1 figures.
        The residual that the toroidal coupling leaves (middle column) has no
        dominant low-degree part (`tables/v2v1.md`).
        """,
        source = @__FILE__,
        runs = [tr.path, nt.path, v.path])
end

# =============================================================================
# Drive
# =============================================================================
v2v1 = load_v2v1()

# Every block D run, for the shared colour scale -- both pairs and the
# lateral-signal figure have to land on one range. See FigStyle.deglac_scales.
scale_runs = Dict(id => load_run(id) for id in SCALE_RUNS)
filter!(kv -> kv.second !== nothing, scale_runs)
SC = deglac_scales(scale_runs, rsl_anomaly, uplift_rate)
@info @sprintf("shared colour limits: rsl ±%.0f m, Δrsl ±%.0f m, uplift ±%.1f, Δuplift ±%.2f mm/yr",
               SC.rsl[2], SC.drsl[2], SC.upl[2], SC.dupl[2])

series = Pair{String,Any}[]
open(joinpath(tables_dir(), "v2v1.md"), "w") do io
    println(io, "# VILMA2 against VILMA1 — the last deglaciation, ℓmax 170\n")
    println(io, "Generated by `analysis/check_v2v1.jl`.\n")
    println(io, "Both backends run behind the same driver on the same GLAC-1D forcing ",
                "and the same output grid, referenced to the same file — VILMA2's ",
                "own LGM reference state (`i_eq = 2`; the 26 ka load held 30 kyr) — with the degree-1 frame ",
                "matched (`deg1_frame = \"cm\"`). The viscosity field is shared by ",
                "construction: `experiments/make_vilma1_visc.jl` applies VILMA2's own ",
                "clamp to the Bagge et al. (2021) file both backends read ",
                "(`validation/data/input/bagge2021.nc`, lon-lat-r, not the FesmData ",
                "product that is the model default), because VILMA1 reads its ",
                "NetCDF directly through `io_nc3in` and none of the native clamping reaches ",
                "it. `C_ocean` and `bsl` are VILMA2 diagnostics computed the same way ",
                "for both backends and are therefore not independent between them.\n")
    println(io, "**With an LGM datum there is no initialization asymmetry to caveat.** ",
                "VILMA1 begins at its own reference with no first-step jump and VILMA2 ",
                "begins at a held-load state of its own Earth model. Neither needs a spin-up. ",
                "The `i_eq = 1` pairs, and the seven datum-hunting pairs that preceded ",
                "them, are retired.\n")

    for p in PAIRS
        f = load_run(p.v2)
        v = load_run(p.v1)
        if f === nothing || v === nothing
            @warn "pair $(p.tag) is not published -- skipping" p.v2 p.v1
            continue
        end
        @info "pair $(p.tag): V2 $(length(f.time)) steps, V1 $(length(v.time)) steps"

        maps_figure(p, f, v, SC)

        # The residual decomposition behind the summary figure. Not recomputed
        # here: every one of its 27 epochs needs both runs' full fields
        # differenced and a degree-0/1 fit, which the published epoch subset
        # cannot supply, so make_reduced.jl computes it once and publishes it.
        k = findfirst(x -> x.tag == p.tag, v2v1)
        if k === nothing
            @warn "pair $(p.tag) has no published residual series -- skipping the summary"
        else
            push!(series, p.tag => v2v1[k])
        end

        # ------------------------------------------------------------ table
        println(io, "\n## Pair `$(p.tag)` — $(p.lab)\n")
        println(io, "| time | rms [m] | max\\|Δ\\| [m] | 99th pct [m] | V2 signal max [m] | rms/signal | deg 0 [m] | \\|deg 1\\| [m] | rms after deg 0+1 |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for t in T_MARKS
            d = rsl_anomaly(f, t) .- rsl_anomaly(v, t)
            sg = maximum(abs.(rsl_anomaly(f, t))); st = stats(d)
            c0, c1, r, _ = deg01(d, f.lon, f.lat)
            @printf(io, "| %.0f ka | %.3f | %.3f | %.3f | %.1f | %.2e | %.2f | %.2f | %.3f |\n",
                    ka(t), st.rms, st.mx, st.p99, sg, st.rms/sg, c0, c1, sqrt(mean(r .^ 2)))
        end

        dupl = uplift_rate(f) .- uplift_rate(v)
        st = stats(dupl); sg = maximum(abs.(uplift_rate(f)))
        @printf(io, "\n**Present-day uplift rate**: rms **%.4f mm/yr**, max **%.4f**, 99th pct **%.4f**, ",
                st.rms, st.mx, st.p99)
        @printf(io, "against a V2 signal max of %.2f mm/yr (rms/signal %.2e).\n", sg, st.rms/sg)

        bf = Float64.(f.ds["bsl"][:]); bv = Float64.(v.ds["bsl"][:])
        nb = min(length(bf), length(bv))
        @printf(io, "\n**Barystatic sea level**: max \\|V2 − V1\\| = **%.3f m**; at present day %.3f m.\n",
                maximum(abs.(bf[1:nb] .- bv[1:nb])), bf[nb] - bv[nb])
        omf, omv = ocean_mean(f, 0.0), ocean_mean(v, 0.0)
        @printf(io, "\n**Ocean-mean rsl at present day**: V2 %.2f m, V1 %.2f m, difference %.2f m.\n",
                omf, omv, omf - omv)
        if k !== nothing
            sr = v2v1[k]; ok = findall(isfinite, sr.rms); i = ok[argmax(sr.rms[ok])]
            okr = findall(isfinite, sr.rel); j = okr[argmax(sr.rel[okr])]
            @printf(io, "\n**Peak of the residual series** (`deglac-v2v1-summary.png`, every 1 kyr): rms %.2f m at %.0f ka; rms/signal at most %.1f %% (%.0f ka).\n",
                    sr.rms[i], -sr.time[i] / 1000, 100 * sr.rel[j], -sr.time[j] / 1000)
            record!("v2v1", "$(p.tag)_rms_peak_m" => round(sr.rms[i], digits = 2),
                            "$(p.tag)_rms_peak_ka" => round(-sr.time[i] / 1000, digits = 1))
        end

        println(io, "\n| time | ocean-function cells differing | % of grid |")
        println(io, "|---|---:|---:|")
        for t in [T_MARKS; 0.0]
            c1 = slice_at(f, "C_ocean", t); c2 = slice_at(v, "C_ocean", t)
            nd = count(abs.(c1 .- c2) .> 0.5)
            @printf(io, "| %.0f ka | %d | %.3f %% |\n", ka(t), nd, 100nd/length(c1))
        end

        d21 = rsl_anomaly(f, -21000.0) .- rsl_anomaly(v, -21000.0)
        c0_21, c1_21, _, _ = deg01(d21, f.lon, f.lat)
        record!("v2v1",
            "$(p.tag)_deg0_21ka_m"     => round(c0_21, digits = 2),
            "$(p.tag)_deg1_21ka_m"     => round(c1_21, digits = 2),
            "$(p.tag)_ocean_mean_F_m"  => round(omf, digits = 2),
            "$(p.tag)_ocean_mean_V_m"  => round(omv, digits = 2),
            "$(p.tag)_ocean_mean_dif_m" => round(omf - omv, digits = 2),
            "$(p.tag)_rms_21ka_m"      => round(sqrt(mean(d21 .^ 2)), digits = 2),
            "$(p.tag)_rel_21ka_pct"    => round(100 * sqrt(mean(d21 .^ 2)) / maximum(abs.(rsl_anomaly(f, -21000.0))), digits = 1),
            "$(p.tag)_uplift_rms_mmyr" => round(st.rms, digits = 3),
            "$(p.tag)_bsl_max_dif_m"   => round(maximum(abs.(bf[1:nb] .- bv[1:nb])), digits = 3))
    end

    # ------------------------------------------------- degree 2, order 1
    println(io, "\n## The degree-2, order-1 part of the residual at 21 ka\n")
    println(io, "Area-weighted (Gauss quadrature) rms of V2 − V1 relative sea level at 21 ka, ",
                "its Y21 part (cos and sin, by projection), the rms after removing it, and ",
                "the largest single (degree, order) part up to degree 8 after Y21. At ",
                "907bc74 the table read 3-D 8.91 / 7.74 / 4.41 m and 1Db 7.93 / 7.86 / 1.04 m.\n")
    println(io, "| pair | rms [m] | Y21 part [m] | rms without Y21 [m] | largest other part [m] |")
    println(io, "|---|---:|---:|---:|---|")
    for p in PAIRS[1:2]
        f, v = load_run(p.v2), load_run(p.v1)
        (f === nothing || v === nothing) && continue
        d = rsl_anomaly(f, -21000.0) .- rsl_anomaly(v, -21000.0)
        tot, parts = sh_parts(d, f.lon, f.lat)
        y21, fit21 = parts[(2, 1)]
        W = repeat(reshape(lat_weights(f.lat), 1, :), length(f.lon), 1) ./ length(f.lon)
        rest = sqrt(sum(W .* (d .- fit21) .^ 2))
        others = sort([(k, v[1]) for (k, v) in parts if k != (2, 1)], by = x -> -x[2])
        (lk, mk), rk = others[1]
        @printf(io, "| %s | %.2f | %.2f | %.2f | degree %d, order %d: %.2f |\n",
                p.tag, tot, y21, rest, lk, mk, rk)
        record!("v2v1", "$(p.tag)_rms_area_21ka_m" => round(tot, digits = 2),
                        "$(p.tag)_y21_21ka_m" => round(y21, digits = 2),
                        "$(p.tag)_rms_noy21_21ka_m" => round(rest, digits = 2),
                        "$(p.tag)_largest_other_lm" => "$(lk),$(mk)",
                        "$(p.tag)_largest_other_m" => round(rk, digits = 2))
    end
    println(io, "\nParts are projections on the real harmonics of the grid; the rms of ",
                "everything above degree 8 is what the first and last columns leave.")

    println(io, "\n## Differences that are not solver error\n")
    println(io, "* **Load path within a coupling interval** — VILMA2 ramps the load ",
                "linearly across the 100 yr interval and sub-steps adaptively to its own ",
                "Maxwell stability limit; VILMA1 holds the load at the interval value and ",
                "takes `vilma_v1_nsub` uniform sub-steps (12 here, an 8.33 yr step against a ",
                "measured 11.0 yr Maxwell limit for this structure).")
    println(io, "* **Radial discretisation** — the two codes integrate the same viscosity ",
                "field onto their own element grids.")
    println(io, "* **Degree-1 frame** — matched here (`cm` in both), and panel (c) of the ",
                "summary figure is what keeps it checkable.")

    # ---------------------------------------------------------------- V7
    tr, nt, v = load_run(TOR), load_run(NOTOR), load_run(V3D)
    if any(isnothing, (tr, nt, v))
        @warn "V7 runs not published -- skipping the toroidal section" TOR NOTOR V3D
    else
        println(io, "\n## V7 — the toroidal coupling\n")
        println(io, "`$(TOR)` and `$(NOTOR)`, one binary, with and without the ",
                    "toroidal coupling; V2−V1 against `$(V3D)`. rms over the grid.\n")
        println(io, "| time | tor − notor rms [m] | max\\|Δ\\| [m] | V2−V1 notor [m] | V2−V1 tor [m] |")
        println(io, "|---|---:|---:|---:|---:|")
        for t in T_MARKS
            d = rsl_anomaly(tr, t) .- rsl_anomaly(nt, t)
            @printf(io, "| %.0f ka | %.3f | %.2f | %.3f | %.3f |\n", ka(t),
                    stats(d).rms, stats(d).mx,
                    stats(rsl_anomaly(nt, t) .- rsl_anomaly(v, t)).rms,
                    stats(rsl_anomaly(tr, t) .- rsl_anomaly(v, t)).rms)
        end
        du  = uplift_rate(tr) .- uplift_rate(nt)
        d_notor = stats(uplift_rate(nt) .- uplift_rate(v)).rms
        d_tor = stats(uplift_rate(tr) .- uplift_rate(v)).rms
        @printf(io, "\n**Present-day uplift rate**: tor − notor rms **%.4f mm/yr**, max %.3f; ",
                stats(du).rms, stats(du).mx)
        @printf(io, "V2−V1 rms %.4f (notor) → %.4f mm/yr (tor).\n", d_notor, d_tor)

        d21 = rsl_anomaly(tr, -21000.0) .- rsl_anomaly(nt, -21000.0)
        record!("v7",
            "tor_minus_notor_rms_21ka_m"  => round(stats(d21).rms, digits = 2),
            "tor_minus_notor_uplift_rms_mmyr" => round(stats(du).rms, digits = 3),
            "v2v1_uplift_rms_notor_mmyr"    => round(d_notor, digits = 3),
            "v2v1_uplift_rms_tor_mmyr"      => round(d_tor, digits = 3))
        toroidal_figure(tr, nt, v, SC)
    end
end
@info "wrote validation/tables/v2v1.md"

isempty(series) || summary_figure(series)
