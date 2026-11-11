#!/usr/bin/env julia
#
# v2v1_series_maps.jl — the VILMA2–VILMA1 deglaciation comparison as pictures:
# rsl through the deglaciation as time series, and as maps with their difference.
#
#   julia --project=. analysis/v2v1_series_maps.jl
#
# WHAT IS PLOTTED, and why it is the RAW rsl rather than the present-day-
# referenced anomaly. `rsl` is already measured against each run's own reference
# state, which is the LGM equilibrium for the i_eq=2 pairs and RTopo-2 for the
# i_eq=1 pair (stated per panel). Subtracting the present-day field on top of
# that — as the maps in check_v2v1.jl do — hides the shape of the disagreement:
# it makes the V2−V1 residual look as though it DECAYS from 94 m at 21 ka to 5 m
# at 2 ka, when in the raw fields it GROWS monotonically from 3.6 m to 88 m.
#.
#
# The summary curves are chosen to separate the two things that differ:
#   ocean-mean rsl  the mass/eustatic term -- where the ~60 m offset lives in
#                   the i_eq=2 pairs, and where the two codes agree to 0.6 m
#                   in the i_eq=1 pair
#   pattern rms     each field about its OWN ocean mean: the fingerprint, whose
#                   amplitude differs ~3x and whose correlation goes NEGATIVE
#                   late in the i_eq=2 pairs
#   max|rsl|        the headline scalar the run log prints

using CairoMakie
using Printf
using Statistics

include(joinpath(@__DIR__, "Published.jl"))
include(joinpath(@__DIR__, "FigStyle.jl"))
using .Published
using .FigStyle

const REPO = dirname(@__DIR__)

# tag, V2 dir, V1 dir, datum, note
# The two pairs of the pruned list: both on
# VILMA2's own LGM equilibrium (i_eq = 2) with the degree-1 frames matched
# (deg1_frame = "cm"). See analysis/check_v2v1.jl for why those two choices, and
# why the seven datum-hunting pairs this table used to carry are gone.
const PAIRS = [
    ("3D",  "3d",  "D-V2-3D.l170.fe",  "D-V1-3D.l170.fe",
     "LGM equilibrium (i_eq=2)", "3-D Bagge viscosity; degree-1 frames matched"),
    ("1Db", "1db", "D-V2-1Db.l170.fe", "D-V1-1Db.l170.fe",
     "LGM equilibrium (i_eq=2)", "radial mean of the same field; isolates the lateral signal"),
]

const EPOCHS = [-21000.0, 0.0]        # maps are drawn at these
const C_D = :black

# The three curves and the V2−V1 residual are read, not recomputed: each one needs
# every time level of an 850 MB file, which is exactly the kind of derivation the
# published epoch subset cannot serve, so make_reduced.jl computes them once. The
# area weighting that goes into (b) and (c) is cos(lat) and lives there; see
# Published.raw_series for why it is deliberately not the Gaussian quadrature
# weighting used by `ocean_mean`.

# ---------------------------------------------------------------- time series
function do_series(tag, key, fid, vid, datum, note, v2v1)
    f = load_run(fid); v = load_run(vid)
    (f === nothing || v === nothing) && (@warn "series $tag not published -- skipping"; return false)
    k = findfirst(x -> x.tag == key, v2v1)
    k === nothing && (@warn "no published V2−V1 residual for $tag -- skipping"; return false)
    idx = findall(t -> t >= -26000.0, f.time)[1:2:end]    # every 200 yr is plenty
    ts = f.time[idx]
    sf = Published.raw_series(f); sv = Published.raw_series(v)
    sf = (mx = sf.max_abs[idx], om = sf.omean[idx], pr = sf.pattern[idx])
    sv = (mx = sv.max_abs[idx], om = sv.omean[idx], pr = sv.pattern[idx])
    dr = v2v1[k].raw_rms[idx]
    kyr = ts ./ 1000

    fig = Figure(size = (1100, 640))

    ax1 = Axis(fig[1, 1], ylabel = "max |rsl| [m]", title = "peak signal")
    lines!(ax1, kyr, sf.mx; color = C_V2, label = "VILMA2")
    lines!(ax1, kyr, sv.mx; color = C_V1, label = "VILMA1")
    axislegend(ax1, position = :lt)

    ax2 = Axis(fig[1, 2], ylabel = "ocean-mean rsl [m]", title = "mass / eustatic term")
    lines!(ax2, kyr, sf.om; color = C_V2); lines!(ax2, kyr, sv.om; color = C_V1)

    ax3 = Axis(fig[2, 1], xlabel = "time [ka]", ylabel = "rms about own ocean mean [m]",
               title = "fingerprint pattern")
    lines!(ax3, kyr, sf.pr; color = C_V2); lines!(ax3, kyr, sv.pr; color = C_V1)

    ax4 = Axis(fig[2, 2], xlabel = "time [ka]", ylabel = "V2 − V1 [m]", title = "disagreement")
    lines!(ax4, kyr, dr; color = C_D, label = "rms of V2 − V1")
    lines!(ax4, kyr, sf.om .- sv.om; color = C_GREY, label = "ocean-mean V2 − V1")
    hlines!(ax4, [0.0]; color = (:black, 0.25), linestyle = :dash, linewidth = 0.8)
    axislegend(ax4, position = :lb)

    panels!(fig, [(1, 1), (1, 2), (2, 1), (2, 2)])
    rowgap!(fig.layout, 12); colgap!(fig.layout, 18)

    savefig(fig, "v2v1-raw-series-$(lowercase(tag))";
        title = "V2 − V1 on the raw sea-level field — $tag",
        caption = """
        VILMA2 (red) and VILMA1 (blue) through the last deglaciation for the $tag Earth
        structure, on the RAW relative-sea-level field rather than the present-day anomaly.
        (a) Peak |rsl| over the globe. (b) Relative sea level averaged over the ocean of
        the moment. (c) Root-mean-square of the field about its own ocean mean — the
        fingerprint, with the mass term removed. (d) The disagreement: rms of V2 − V1 over
        the grid, and the difference of the two ocean means. rsl is measured against
        $datum.
        """,
        notes = """
        Diagnostic. This is the RAW field, which is why it is kept separate from the paper's
        V2−V1 figures: those subtract each run's present-day field, and doing so reverses the
        apparent time evolution of the residual. Against the present-day datum the residual
        looks as though it DECAYS through the deglaciation; on the raw field it GROWS
        monotonically, because the datum itself carries most of the disagreement. Both are
        true statements about different quantities, and confusing them cost a session.
        Panel (c) is the diagnostic that separates a mass disagreement from a pattern
        disagreement — the two codes can move the same amount of water into a different
        shape, and (b) alone cannot see that.
        """,
        source = @__FILE__,
        runs = [f.path, v.path, Published.v2v1_path()])
    true
end

# ---------------------------------------------------------------- maps
function do_maps(tag, fid, vid, datum, note)
    f = load_run(fid); v = load_run(vid)
    (f === nothing || v === nothing) && (@warn "maps $tag not published -- skipping"; return false)
    for t in EPOCHS
        a = slice_at(f, "rsl", t); b = slice_at(v, "rsl", t); d = a .- b
        lim  = (q = quantile(abs.(vec([a; b])), 0.995); (-q, q))
        dlim = (q = quantile(abs.(vec(d)), 0.995); (-q, q))
        fig = Figure(size = (1450, 420))
        for (i, (fld, lab, cr)) in enumerate(((a, "VILMA2", lim), (b, "VILMA1", lim),
                                              (d, "VILMA2 − VILMA1", dlim)))
            ax = Axis(fig[1, i]; title = lab, aspect = DataAspect(),
                      xticks = ([-90, 0, 90], ["90°W", "0°", "90°E"]),
                      yticks = ([-60, 0, 60], ["60°S", "0°", "60°N"]),
                      yticklabelsvisible = i == 1)
            limits!(ax, -180, 180, -90, 90)
            hm = heatmap!(ax, f.lon, f.lat, fld; colormap = CMAP_DIV, colorrange = cr)
            Colorbar(fig[2, i], hm, vertical = false, flipaxis = false, height = 11,
                     label = "rsl [m]")
        end
        panels!(fig, [(1, 1), (1, 2), (1, 3)])
        rowgap!(fig.layout, 6); colgap!(fig.layout, 8)

        ka = round(Int, -t/1000)
        savefig(fig, @sprintf("v2v1-raw-map-%s-%dka", lowercase(tag), ka);
            title = "Raw sea-level field at $ka ka — $tag",
            caption = """
            The raw relative-sea-level field at $ka ka for the $tag Earth structure, measured
            against $datum rather than against each run's present day. (a) VILMA2.
            (b) VILMA1, on the same colour scale. (c) Their difference, on its own scale.
            """,
            notes = """
            Diagnostic, and the companion to the raw time series. The paper's V2−V1 map figures
            show the present-day anomaly, which is the right quantity for comparing against
            observations but hides where the disagreement between the two codes actually
            sits: most of it is in the datum. These panels show the fields the solvers
            produced.
            """,
            source = @__FILE__,
            runs = [f.path, v.path])
    end
    true
end

let done = 0, v2v1 = load_v2v1()
    for (tag, key, fid, vid, datum, note) in PAIRS
        do_series(tag, key, fid, vid, datum, note, v2v1) &&
            do_maps(tag, fid, vid, datum, note) && (done += 1)
    end
    @info "completed $done of $(length(PAIRS)) pairs"
end
