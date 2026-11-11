#!/usr/bin/env julia
#
# run_blockD_diags.jl — the block D diagnostics that are NOT maps of one run
#.
#
#   julia --project=. analysis/run_blockD_diags.jl
#
# Reads `data/published/blockD/`, not `runs/`. The sea-level budget is not
# recomputed here — it needs every time level, so it is computed once by
# analysis/make_reduced.jl and published as a series.
#
# Produces three numbered figures and the block D table:
#
#   barystatic-and-rotation  the global scalars: barystatic sea level, the
#                            sea-level budget against the ice actually lost, and
#                            the polar motion the rotational feedback produces
#   lateral-signal           3-D minus 1Db in maps and in present-day uplift —
#                            the figure the paper's lateral-viscosity argument
#                            rests on
#   site-time-series         relative sea level at eight sites spanning
#                            near-field, Antarctic, far-field and intermediate
#                            settings
#   validation/tables/blockD.md   the difference table
#
# WHY 3-D MINUS 1Db AND NOT 3-D MINUS 1Da: 1Db is the radial mean OF the same
# Bagge field, so the difference isolates the LATERAL variation. 1Da would be a
# different profile entirely, so 3-D minus 1Da conflates lateral structure with
# a change in the radial reference. The 1Da row is not in the pruned list.

using CairoMakie
using NCDatasets
using Printf
using Statistics

include(joinpath(@__DIR__, "Published.jl"))
include(joinpath(@__DIR__, "FigStyle.jl"))
using .Published
using .FigStyle

const REPO = dirname(@__DIR__)

# The production matrix: two Earth structures, two solvers, one truncation.
# Structure is carried by COLOUR and solver by LINE STYLE, so neither axis needs
# a sentence of legend to be read off the figure.
#
# `lw` is not decoration. On every global scalar here the four runs agree to
# within the thickness of a line — which is the result — so the 3-D curves are
# drawn wide and translucent and the 1Db curves thin on top of them. Without
# that the figure shows one curve and silently asserts that the others were
# never run.
const SERIES = [
    (id = "D-V2-3D.l170.fe",  lab = "VILMA2, 3-D", col = C_3D, ls = :solid, lw = 3.4, al = 0.5),
    (id = "D-V2-1Db.l170.fe", lab = "VILMA2, 1Db", col = C_1D, ls = :solid, lw = 1.4, al = 1.0),
    (id = "D-V1-3D.l170.fe",  lab = "VILMA1, 3-D",       col = C_3D, ls = :dash,  lw = 3.4, al = 0.5),
    (id = "D-V1-1Db.l170.fe", lab = "VILMA1, 1Db",       col = C_1D, ls = :dash,  lw = 1.4, al = 1.0),
]

const F3D, F1D = "D-V2-3D.l170.fe", "D-V2-1Db.l170.fe"

# The eight sites, from the published schema rather than from a list here: the
# series are extracted during reduction, so asking for a site that was not
# reduced is a request to re-run make_reduced.jl, and Published.site_series says
# so rather than returning the wrong cell.
const SITES = Published.SITES

const T_MARKS = [-21000.0, -14000.0, -8000.0, 0.0]

# ---------------------------------------------------------------- load
# No completeness check: make_reduced.jl refuses to publish a run that has not
# reached the present day, so a file that exists is a finished run.
runs = Dict{String,Any}()
for s in SERIES
    r = load_run(s.id)
    r === nothing ? (@warn "$(s.id) is not in data/published/blockD — skipping") :
                    (runs[s.id] = r)
end
isempty(runs) && error("no published block D runs; run ./run.sh reduce")
@info "loaded $(length(runs)) run(s): $(join(sort(collect(keys(runs))), ", "))"
avail = [s for s in SERIES if haskey(runs, s.id)]
have_pair = haskey(runs, F3D) && haskey(runs, F1D)

# =============================================================================
# Fig — barystatic sea level, budget closure, polar motion
# =============================================================================
fig = Figure(size = (1150, 720))

ax1 = Axis(fig[1, 1], xlabel = "time [ka]", ylabel = "barystatic sea level [m]")
for s in avail
    r = runs[s.id]
    lines!(ax1, r.time ./ 1000, Float64.(r.ds["bsl"][:]);
           color = (s.col, s.al), linestyle = s.ls, linewidth = s.lw, label = s.lab)
end
axislegend(ax1, position = :lt)

# The differences, because panel (a) cannot show them: the four curves agree to
# a fraction of a metre on a 105 m signal, so on (a)'s axis they are one line.
ax2 = Axis(fig[1, 2], xlabel = "time [ka]", ylabel = "Δ barystatic sea level [m]")
bsl(id) = Float64.(runs[id].ds["bsl"][:])
tvec(id) = runs[id].time ./ 1000
for (a, b, lab, c, ls) in (("D-V2-3D.l170.fe",  "D-V1-3D.l170.fe",  "V2 − V1, 3-D",  C_3D, :solid),
                           ("D-V2-1Db.l170.fe", "D-V1-1Db.l170.fe", "V2 − V1, 1Db",  C_1D, :solid),
                           ("D-V2-3D.l170.fe",  "D-V2-1Db.l170.fe", "3-D − 1Db, V2", :black, :dot))
    (haskey(runs, a) && haskey(runs, b)) || continue
    x, y = bsl(a), bsl(b); n = min(length(x), length(y))
    lines!(ax2, tvec(a)[1:n], x[1:n] .- y[1:n]; color = c, linestyle = ls, label = lab)
end
hlines!(ax2, [0.0], color = (:black, 0.25), linestyle = :dash, linewidth = 0.8)
axislegend(ax2, position = :lb)

# Budget closure: what the model delivered against what the ice released. The
# target is computed from the run's own reference ice and bed, so the two curves
# are commensurate by construction and the gap between them is the claim.
ax3 = Axis(fig[2, 1], xlabel = "time [ka]", ylabel = "ocean-mean rsl [m]")
target_drawn = false
for s in avail
    startswith(s.id, "D-V2") || continue
    b = budget_series(runs[s.id])
    lines!(ax3, b.time ./ 1000, b.rsl_mean;
           color = (s.col, s.al), linewidth = s.lw, label = s.lab * " delivered")
    if !target_drawn
        lines!(ax3, b.time ./ 1000, b.target; color = :black, linestyle = :dot,
               linewidth = 1.4, label = "above-flotation volume lost")
        global target_drawn = true
    end
end
axislegend(ax3, position = :lt)

# Polar motion. VILMA2 only: the VILMA1 backend does not write the rotation
# diagnostic, so drawing an empty V1 line would invite the reader to conclude
# the two codes disagree about polar motion when they were never compared.
ax4 = Axis(fig[2, 2], xlabel = "time [ka]", ylabel = "polar motion [10⁻³]")
for s in avail
    startswith(s.id, "D-V2") || continue
    pm = polar_motion(runs[s.id])
    pm === nothing && continue
    m1, m2 = pm
    t = runs[s.id].time ./ 1000
    lines!(ax4, t, m1 .* 1e3; color = (s.col, s.al), linewidth = s.lw, label = s.lab * ", m₁")
    lines!(ax4, t, m2 .* 1e3; color = (s.col, s.al), linewidth = s.lw,
           linestyle = :dashdot, label = s.lab * ", m₂")
end
axislegend(ax4, position = :lb)

panels!(fig, [(1, 1), (1, 2), (2, 1), (2, 2)])
colgap!(fig.layout, 18); rowgap!(fig.layout, 14)

savefig(fig, "barystatic-and-rotation";
    title = "Global scalars of the deglaciation",
    caption = """
    Globally integrated diagnostics of the last deglaciation at ℓmax 170, for both
    Earth structures and both solid-Earth solvers. Colour distinguishes the Earth
    structure (red, the 3-D Bagge et al. (2021) field; blue, 1Db, its radial mean);
    line style distinguishes the solver (solid, VILMA2; dashed, VILMA1); the
    3-D curves are drawn wide and translucent so that the 1Db curves, which lie on
    top of them, remain visible. (a) Barystatic sea level, the model's own eustatic
    bookkeeping, rising from the LGM datum. (b) The same quantity as differences:
    VILMA2 minus VILMA1 for each structure, and 3-D minus 1Db within VILMA2.
    (c) Sea-level budget closure for VILMA2: delivered ocean-mean relative sea
    level against the above-flotation volume of the reference ice lost by that time
    (black dotted, shown for the 3-D run), both in barystatic metres. (d) Polar
    motion components m₁ and m₂ (dash-dotted) from the rotational feedback,
    VILMA2 only.
    """,
    notes = """
    The four curves in (a) agree to well under a metre on a ~105 m signal, which is
    why (b) exists: on (a)'s axis the result "the two solvers and the two structures
    move the same amount of water" is indistinguishable from "only one run was made".
    Panel (c) is the check that caught the marine-reference mask defect (model
    `5385a8b`): the target curve is rebuilt from the run's own reference
    ice and bed and is never seen by the sea-level equation, so it can contradict
    what the solver produced. Panel (d) carries no VILMA1 line because the VILMA1
    backend does not write the rotation diagnostic — the two codes are not being
    compared here.
    Panels (a) and (c) measure different things and should not be read against each
    other: `bsl` is eustatic bookkeeping, the (c) curve is an area average of the
    solved sea-level field over the ocean of the moment.
    """,
    source = @__FILE__,
    runs = [runs[s.id].path for s in avail])

# =============================================================================
# Fig — the lateral-viscosity signal: 3-D minus 1Db
# =============================================================================
if have_pair
    r3, r1 = runs[F3D], runs[F1D]
    dmaps = [rsl_anomaly(r3, t) .- rsl_anomaly(r1, t) for t in T_MARKS[1:3]]
    dupl  = uplift_rate(r3) .- uplift_rate(r1)
    # Both scales are needed, for two different questions. The SHARED range
    # (FigStyle.deglac_scales, computed over the whole block D matrix) answers
    # "how big is this against everything else"; its OWN range answers "what
    # shape is it". The autoscaled one is the paper figure, the shared one the
    # diagnostic, and the measured ratio to the VILMA2−VILMA1 residual is
    # stated in the notes in whichever direction it falls.
    SC = deglac_scales(runs, rsl_anomaly, uplift_rate)
    dl, ul = SC.drsl, SC.dupl

    # Drawn twice, from one function: autoscaled as the paper figure, and on the
    # shared scale as a diagnostic. Two figures rather than one because they
    # answer different questions — what shape is this, and how big is it against
    # everything else — and a single autoscaled figure with no companion silently
    # answers the second one wrongly, which is what the set used to do.
    function lateral_figure(slug, dl, ul, scale_note, companion, levels)
        fig3 = Figure(size = (1200, 620))
        hm = hu = nothing
        for (k, (t, dm)) in enumerate(zip(T_MARKS[1:3], dmaps))
            row, col = divrem(k - 1, 2) .+ (1, 1)
            a = Axis(fig3[row, col], title = @sprintf("%.0f ka", ka(t)),
                     aspect = DataAspect(),
                     xticks = ([-90, 0, 90], ["90°W", "0°", "90°E"]),
                     yticks = ([-60, 0, 60], ["60°S", "0°", "60°N"]),
                     xticklabelsize = TICK_FS, yticklabelsize = TICK_FS,
                     xticklabelsvisible = row == 2, yticklabelsvisible = col == 1)
            limits!(a, -180, 180, -90, 90)
            hm = heatmap!(a, r3.lon, r3.lat, dm; colormap = CMAP_DIV, colorrange = dl)
            contour!(a, r3.lon, r3.lat, dm; levels, color = (:black, 0.45), linewidth = 0.5)
        end
        a4 = Axis(fig3[2, 2], title = "present day, uplift rate", aspect = DataAspect(),
                  xticks = ([-90, 0, 90], ["90°W", "0°", "90°E"]),
                  yticks = ([-60, 0, 60], ["60°S", "0°", "60°N"]),
                  xticklabelsize = TICK_FS, yticklabelsize = TICK_FS,
                  yticklabelsvisible = false)
        limits!(a4, -180, 180, -90, 90)
        hu = heatmap!(a4, r3.lon, r3.lat, dupl; colormap = :curl, colorrange = ul)
        Colorbar(fig3[3, 1], hm, vertical = false, flipaxis = false,
                 label = "Δ relative sea level [m]", height = 11)
        Colorbar(fig3[3, 2], hu, vertical = false, flipaxis = false,
                 label = "Δ uplift rate [mm/yr]", height = 11)
        panels!(fig3, [(1, 1), (1, 2), (2, 1), (2, 2)])
        rowgap!(fig3.layout, 6); colgap!(fig3.layout, 10)

        savefig(fig3, slug;
            title = "The lateral-viscosity signal",
            caption = """
            The difference made by lateral viscosity structure: the 3-D Bagge et al. (2021)
            run minus the 1Db run, which is the radial mean of that same field, at ℓmax 170.
            (a–c) Relative sea level relative to each run's own present day, at 21, 14 and
            8 ka, with contours at $(join(filter(>(0), levels), " and ")) m and their
            negatives. (d) Present-day bedrock uplift rate. $scale_note
            """,
            notes = """
            1Db is the radial mean of the same field rather than an independent 1-D profile,
            so this difference isolates lateral variation instead of mixing it with a change
            of radial reference. There is no 0 ka sea-level panel on purpose: the anomaly is
            measured against each run's own present day, so at 0 ka it is identically zero
            for both and the difference would carry no information. Present day is covered
            by the uplift rate in (d).
            THE LATERAL SIGNAL AGAINST THE SOLVER RESIDUAL. Over 21–8 ka it is
            $(@sprintf("%.2f", lat_rms)) m rms in relative sea level, against
            $(@sprintf("%.1f", v2v1_rms)) m rms for VILMA2 minus VILMA1 on the same
            epochs: $(lat_rms >= v2v1_rms ? "larger" : "smaller") by a factor of
            $(@sprintf("%.1f", max(lat_rms, v2v1_rms) / min(lat_rms, v2v1_rms))). An
            autoscaled colour bar fills its range whatever the magnitude, so that
            comparison is read from the numbers, or by eye from the companion figure
            `$companion`, which is this same difference on the other scale.
            """,
            source = @__FILE__,
            runs = [runs[F3D].path, runs[F1D].path])
    end

    # The two rms values the caption weighs against each other, over the epochs
    # actually drawn, so the claim in the notes is measured rather than recalled.
    rmsof(fs) = sqrt(mean(vcat([vec(f) for f in fs]...) .^ 2))
    lat_rms = rmsof(dmaps)
    v2v1_rms = let acc = Matrix{Float64}[]
        for (a, b) in (("D-V2-3D.l170.fe", "D-V1-3D.l170.fe"),
                       ("D-V2-1Db.l170.fe", "D-V1-1Db.l170.fe"))
            (haskey(runs, a) && haskey(runs, b)) || continue
            append!(acc, [rsl_anomaly(runs[a], t) .- rsl_anomaly(runs[b], t) for t in T_MARKS[1:3]])
        end
        isempty(acc) ? NaN : rmsof(acc)
    end

    lateral_figure("lateral-signal", symlim(dmaps...), symlim(dupl),
        "Colour scales are this figure's own, symmetric about zero and saturating " *
        "at the 99.5th percentile of |difference|, so that the pattern is legible; " *
        "they are therefore NOT comparable with the other block D maps, which share " *
        "one range. Read the magnitude off the scale bar, not off the colours.",
        "figures/diagnostics/lateral-signal-shared-scale",
        [-5, -2, 2, 5])
    lateral_figure("lateral-signal-shared-scale", dl, ul,
        "Colour scales are shared with every other block D map, so a difference " *
        "here and a difference in the VILMA2−VILMA1 figures are directly " *
        "comparable.",
        "figures/Fig03-lateral-signal",
        [-5, -2, 2, 5])
else
    @warn "lateral-signal needs both $F3D and $F1D"
end

# =============================================================================
# Fig — site time series
# =============================================================================
main = [s for s in avail if startswith(s.id, "D-V2")]
if !isempty(main)
    fig4 = Figure(size = (1150, 560))
    for (k, (name, lat, lon)) in enumerate(SITES)
        row, col = divrem(k - 1, 4) .+ (1, 1)
        local a = Axis(fig4[row, col], title = name, titlealign = :center,
                       xlabel = row == 2 ? "time [ka]" : "",
                       ylabel = col == 1 ? "rsl vs present [m]" : "")
        for s in main
            r = runs[s.id]
            v, _, _ = site_series(r, "rsl", lat, lon)
            lines!(a, r.time ./ 1000, v .- v[end]; color = s.col, linestyle = s.ls, label = s.lab)
        end
        hlines!(a, [0.0], color = (:black, 0.25), linestyle = :dash, linewidth = 0.8)
    end
    Legend(fig4[3, 1:4],
           [LineElement(color = s.col, linestyle = s.ls) for s in main],
           [s.lab for s in main], orientation = :horizontal)
    panels!(fig4, [(divrem(k - 1, 4) .+ (1, 1)) for k in 1:length(SITES)])
    rowgap!(fig4.layout, 8)

    savefig(fig4, "site-time-series";
        title = "Relative sea level at eight sites",
        caption = """
        Relative sea level through the last deglaciation at eight locations, for the
        3-D Bagge et al. (2021) viscosity field (red) and 1Db, its radial mean (blue),
        both VILMA2 at ℓmax 170. (a) Hudson Bay and (b) Ångermanälven are
        near-field sites under the Laurentide and Fennoscandian centres of loading;
        (c) Iceland and (f) Patagonia are intermediate-field; (d) Barbados and
        (e) Tahiti are far-field; (g) Amundsen Sea and (h) Ross Sea are Antarctic.
        Each curve is referenced to its own model present day, so all eight panels pass
        through zero at 0 ka. Note the differing vertical scales.
        """,
        notes = """
        Values are the single nearest grid cell, not an average over a neighbourhood: at
        ℓmax 170 a cell is about 0.5°, and near a steep forebulge gradient the choice of
        cell moves the curve by metres. These are model output at a coordinate, not
        predictions of a specific RSL record, and should not be read against
        observational sea-level curves without the usual site-specific care.
        """,
        source = @__FILE__,
        runs = [runs[s.id].path for s in main])
end

# =============================================================================
# The table
# =============================================================================
open(joinpath(tables_dir(), "blockD.md"), "w") do io
    println(io, "# Block D — last deglaciation, structure differences\n")
    println(io, "Generated by `analysis/run_blockD_diags.jl`. GLAC-1D, 26 ka → present ",
                "at 100 yr, `fe` scheme, rotation on, migrating coastlines, ℓmax 170, ",
                "referenced to VILMA2's own LGM reference state (`i_eq = 2`; the 26 ka load held 30 kyr).\n")

    println(io, "## Runs used\n")
    println(io, "| run | structure | ℓmax | steps | bsl at 21 ka [m] | bsl at present [m] |")
    println(io, "|---|---|---:|---:|---:|---:|")
    for s in avail
        r = runs[s.id]
        b = Float64.(r.ds["bsl"][:])
        i21 = nearest_t(r, -21000.0)
        @printf(io, "| `%s` | %s | %d | %d | %.2f | %.2f |\n",
                s.id, r.visc == "3d" ? "3-D" : "1Db", r.lmax, length(r.time), b[i21], b[end])
    end

    println(io, "\n## Sea-level budget closure\n")
    println(io, "Delivered ocean-mean relative sea level at present day against the ",
                "above-flotation volume of the reference ice that was lost, both in ",
                "barystatic metres and both rebuilt from the run's own output. ",
                "`analysis/sle_budget.jl` prints the full decomposition.\n")
    println(io, "| run | delivered [m] | above-flotation target [m] | delivered/target |")
    println(io, "|---|---:|---:|---:|")
    for s in avail
        startswith(s.id, "D-V2") || continue
        b = budget_at(runs[s.id])
        @printf(io, "| `%s` | %.2f | %.2f | %.1f %% |\n", s.id, b.delivered, b.target, 100b.closure)
        record!("blockD", "$(s.id)_delivered_m" => round(b.delivered, digits = 2),
                          "$(s.id)_target_m"    => round(b.target, digits = 2),
                          "$(s.id)_closure_pct" => round(100b.closure, digits = 1))
    end

    if have_pair
        r3, r1 = runs[F3D], runs[F1D]
        println(io, "\n## 3-D minus 1Db — the lateral-viscosity signal\n")
        println(io, "1Db is the radial mean **of the same Bagge (2021) field**, so this ",
                    "difference isolates lateral variation rather than mixing it with a ",
                    "change of radial profile.\n")
        println(io, "The 0 ka row is omitted on purpose: `rsl_anomaly` is measured against ",
                    "each run's own present-day field, so at 0 ka it is identically zero ",
                    "for both and the difference carries no information. Present day is ",
                    "covered by the uplift rate below.\n")
        println(io, "| time | rms [m] | max\\|Δ\\| [m] | 99th pct \\|Δ\\| [m] | signal max\\|rsl\\| [m] | rms/signal |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for t in T_MARKS[1:3]
            d = rsl_anomaly(r3, t) .- rsl_anomaly(r1, t)
            s = maximum(abs.(rsl_anomaly(r3, t)))
            ad = abs.(vec(d))
            @printf(io, "| %.0f ka | %.3f | %.3f | %.3f | %.1f | %.2e |\n",
                    ka(t), sqrt(mean(d .^ 2)), maximum(ad), quantile(ad, 0.99), s,
                    s == 0 ? NaN : sqrt(mean(d .^ 2))/s)
        end
        du = uplift_rate(r3) .- uplift_rate(r1)
        adu = abs.(vec(du))
        @printf(io, "\nPresent-day uplift rate: rms **%.4f mm/yr**, max **%.4f**, 99th pct **%.4f**, ",
                sqrt(mean(du .^ 2)), maximum(adu), quantile(adu, 0.99))
        @printf(io, "against a signal max of %.2f mm/yr.\n", maximum(abs.(uplift_rate(r3))))
        record!("blockD",
                "lateral_uplift_rms_mmyr"    => round(sqrt(mean(du .^ 2)), digits = 3),
                "lateral_uplift_max_mmyr"    => round(maximum(adu), digits = 2),
                "lateral_rsl_rms_21ka_m"     => round(sqrt(mean((rsl_anomaly(r3, -21000.0) .- rsl_anomaly(r1, -21000.0)) .^ 2)), digits = 2))

        b3 = Float64.(r3.ds["bsl"][:]); b1 = Float64.(r1.ds["bsl"][:])
        n = min(length(b3), length(b1))
        @printf(io, "\nBarystatic sea level: max \\|3-D − 1Db\\| = **%.3f m** over the run.\n",
                maximum(abs.(b3[1:n] .- b1[1:n])))

        println(io, "\n### Ocean-function mismatch\n")
        println(io, "| time | cells differing | % of grid |")
        println(io, "|---|---:|---:|")
        for t in T_MARKS
            c3 = slice_at(r3, "C_ocean", t); c1 = slice_at(r1, "C_ocean", t)
            nd = count(abs.(c3 .- c1) .> 0.5)
            @printf(io, "| %.0f ka | %d | %.3f %% |\n", ka(t), nd, 100nd/length(c3))
        end
    end

    println(io, "\n## Polar motion\n")
    for s in avail
        startswith(s.id, "D-V2") || continue
        pm = polar_motion(runs[s.id])
        pm === nothing && (println(io, "* `$(s.id)` — not written by this run."); continue)
        m1, m2 = pm
        @printf(io, "* `%s` — m₁ spans %+.2e … %+.2e, m₂ spans %+.2e … %+.2e (dimensionless).\n",
                s.id, minimum(m1), maximum(m1), minimum(m2), maximum(m2))
    end
    println(io, "\nThe VILMA1 backend does not write the rotation diagnostic, so m₁ and m₂ ",
                "are VILMA2 only and the two codes are not compared on polar motion.")
end
@info "wrote validation/tables/blockD.md"
