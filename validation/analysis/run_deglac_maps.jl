#!/usr/bin/env julia
#
# run_deglac_maps.jl — the last-deglaciation map figure.
#
#   julia --project=. analysis/run_deglac_maps.jl
#
# Reads `data/published/blockD/`, not `runs/`. See analysis/Published.jl — only
# analysis/make_reduced.jl touches the raw output.
#
# Each sea-level panel carries TWO fields in two distinct colour scales:
#
#   base   relative sea level against the model's own present day [m],
#          diverging, with contours — the quantity a GIA model is judged on
#   over   the ice sheet, as surface elevation [km] in its own sequential scale
#          with topographic contours, drawn only where the ice is grounded and
#          thicker than 10 m
#
# so the load and its sea-level response are legible in the same frame without
# one hiding the other.
#
# Equirectangular lon/lat: GeoMakie is not a dependency, and for a global field
# whose features are the ice sheets and their forebulges a plate carrée frame is
# adequate and keeps the analysis environment small.

using CairoMakie
using Printf
using Statistics

include(joinpath(@__DIR__, "Published.jl"))
include(joinpath(@__DIR__, "FigStyle.jl"))
using .Published
using .FigStyle

const REPO = dirname(@__DIR__)
# The block D PRODUCTION matrix, both structures at the production truncation,
# referenced to VILMA2's own LGM equilibrium (i_eq = 2). Two earlier
# defaults were both wrong for this figure and are recorded here so they are not
# re-adopted:
#
#   runs/deglac      equil_time_max = 0, right for timing the transient but it
#                    starts the solid Earth UNRELAXED under LGM ice, so the LGM
#                    panel showed rsl = 0 and far-field sea level ~+1 m instead
#                    of ~-130 m. Deleted.
#   runs/deglac_i1   correctly initialised, but built on the M3-L70-V01
#                    BENCHMARK Earth structure rather than block D's PREM
#                    geometry with the Bagge-derived viscosity. Fine for block P
#                    timing, wrong for a block D physics figure.
#
# Which raw run each published file came from is recorded in its `source_run`
# attribute, so that history stays attached to the data rather than only here.

# Rows, in the order they should appear. The row label names the Earth structure
# and nothing else — the epochs are the columns and the caption carries the rest.
const CASES = [
    ("D-V2-3D.l170.fe",  "3-D"),
    ("D-V2-1Db.l170.fe", "1Db"),
]

const T_LGM = -21000.0
const T_MID = -10000.0

"""
    panel!(fig, pos, run, kind; ...)

One map. `kind` is `:rsl` (with `yr`) or `:uplift`.
"""
function mappanel!(fig, pos, r, kind; yr = 0.0, rsl_lim = nothing, upl_lim = nothing,
                   title = "", show_y = false, show_x = false)
    ax = Axis(fig[pos...]; title, titlealign = :center,
              # No ±180 label: adjacent panels abut, and the 180°E of one lands
              # on top of the 180°W of the next.
              xticks = ([-90, 0, 90], ["90°W", "0°", "90°E"]),
              yticks = ([-60, 0, 60], ["60°S", "0°", "60°N"]),
              xticklabelsize = TICK_FS, yticklabelsize = TICK_FS,
              xticklabelsvisible = show_x, yticklabelsvisible = show_y,
              aspect = DataAspect())
    limits!(ax, -180, 180, -90, 90)

    if kind === :rsl
        base = rsl_anomaly(r, yr)
        hm = heatmap!(ax, r.lon, r.lat, base, colormap = CMAP_DIV, colorrange = rsl_lim)
        # Contours on a log-ish set of levels: the near field is hundreds of
        # metres while the far field is a few, and even levels show only the former.
        contour!(ax, r.lon, r.lat, base,
                 levels = [-400, -200, -100, -50, -20, -5, 5, 20, 50, 100, 200, 400],
                 color = (:black, 0.35), linewidth = 0.4)
        _, ice = srf_elevation(r, yr)
        if any(isfinite, ice)
            # Slightly translucent: under the ice is exactly where the sea-level
            # signal is largest, and an opaque overlay would hide the most
            # informative part of the field.
            heatmap!(ax, r.lon, r.lat, ice ./ 1000, colormap = CMAP_SEQ,
                     colorrange = (0, 4), alpha = 0.72)
            contour!(ax, r.lon, r.lat, replace(ice, NaN => -1e9) ./ 1000,
                     levels = 0.5:0.5:4.0, color = (:black, 0.45), linewidth = 0.35)
        end
        contour!(ax, r.lon, r.lat, slice_at(r, "C_ocean", yr), levels = [0.5],
                 color = (:black, 0.8), linewidth = 0.6)
        return ax, hm
    else
        u  = uplift_rate(r)
        hm = heatmap!(ax, r.lon, r.lat, u, colormap = :curl, colorrange = upl_lim)
        contour!(ax, r.lon, r.lat, u, levels = [-4, -2, -1, 1, 2, 4, 8],
                 color = (:black, 0.35), linewidth = 0.4)
        _, ice = srf_elevation(r, 0.0)
        any(isfinite, ice) && contour!(ax, r.lon, r.lat, replace(ice, NaN => -1e9) ./ 1000,
                                       levels = [0.0], color = (:black, 0.7), linewidth = 0.7)
        contour!(ax, r.lon, r.lat, slice_at(r, "C_ocean", 0.0), levels = [0.5],
                 color = (:black, 0.8), linewidth = 0.6)
        return ax, hm
    end
end

# ---------------------------------------------------------------- load
# Completeness is not rechecked here: `make_reduced.jl` refuses to publish a run
# that has not reached the present day, so anything in data/published/blockD is
# a finished run. A missing file means that case was not reduced.
# Every block D run, for the colour scale; only CASES are drawn. See
# FigStyle.deglac_scales -- loading just the rows drawn would put this figure
# back on a private scale.
scale_runs = Dict(id => load_run(id) for id in SCALE_RUNS)
filter!(kv -> kv.second !== nothing, scale_runs)

runs = Tuple{Any,String}[]
for (id, lab) in CASES
    r = load_run(id)
    if r === nothing
        @warn "$id is not in data/published/blockD — skipping"
        continue
    end
    push!(runs, (r, lab))
end
isempty(runs) && error("no published deglaciation runs; run ./run.sh reduce")
@info "loaded $(length(runs)) deglaciation run(s)"

# Colour limits from FigStyle, computed over the WHOLE block D matrix rather
# than over the two rows drawn here — so a panel in this figure is comparable
# with a panel in the V2−V1 figures, which is the reason the full matrix is
# loaded above even though only two runs appear.
SC = deglac_scales(scale_runs, rsl_anomaly, uplift_rate)
rsl_lim, upl_lim = SC.rsl, SC.upl
@info @sprintf("colour limits (shared): rsl ±%.0f m, uplift ±%.1f mm/yr", rsl_lim[2], upl_lim[2])

# ---------------------------------------------------------------- figure
nrow = length(runs)
fig = Figure(size = (1400, 270 * nrow + 120))
hm_rsl, hm_upl = nothing, nothing
for (i, (r, lab)) in enumerate(runs)
    islast = i == nrow
    _, h1 = mappanel!(fig, (i, 1), r, :rsl; yr = T_LGM, rsl_lim,
                      title = i == 1 ? "21 ka" : "", show_y = true, show_x = islast)
    mappanel!(fig, (i, 2), r, :rsl; yr = T_MID, rsl_lim,
              title = i == 1 ? "10 ka" : "", show_x = islast)
    _, h2 = mappanel!(fig, (i, 3), r, :uplift; upl_lim,
                      title = i == 1 ? "present day, uplift rate" : "", show_x = islast)
    global hm_rsl = h1
    global hm_upl = h2
    Label(fig[i, 0], lab, rotation = π/2, fontsize = 13, tellheight = false)
end
panels!(fig, [(i, j) for i in 1:nrow for j in 1:3])

Colorbar(fig[nrow + 1, 1:2], hm_rsl, vertical = false, flipaxis = false,
         label = "relative sea level vs present day [m]", height = 11)
Colorbar(fig[nrow + 1, 3], hm_upl, vertical = false, flipaxis = false,
         label = "uplift rate [mm/yr]", height = 11)
rowgap!(fig.layout, 6); colgap!(fig.layout, 6)

savefig(fig, "deglaciation-maps";
    title = "Last deglaciation: sea-level response and ice load",
    caption = """
    Relative sea level and present-day bedrock uplift over the last deglaciation,
    forced by GLAC-1D ice from 26 ka to present at a 100 yr coupling interval, at
    ℓmax 170. Top row (a–c): the 3-D Bagge et al. (2021) viscosity field. Bottom row
    (d–f): 1Db, the radial (log-)mean of that same field, which differs from the top
    row only in lateral structure. (a, d) and (b, e) show relative sea level at 21 ka
    and 10 ka relative to each run's own present-day field, with contours at ±5, 20,
    50, 100, 200 and 400 m; the grey overlay is grounded ice-sheet surface elevation,
    contoured every 500 m, drawn where the ice exceeds 10 m. (c, f) show the
    present-day bedrock uplift rate, differenced over the last 500 yr of the run,
    with contours at ±1, 2, 4 and 8 mm/yr. The black line in every panel is the
    coastline implied by the model's own ocean function at that epoch. Colour scales
    are shared across rows.
    """,
    notes = """
    The datum is each run's own present day rather than the reference topography, so
    a non-zero present-day bias cancels and the panels show the change the model
    produced. Both runs are referenced to VILMA2's LGM reference state, the 26 ka
    load held 30 kyr (`i_eq = 2`, `data/reference/vilma2_lgm_eq.nc`), so neither carries a spin-up
    transient. 21 ka rather than the 26 ka start is used throughout block D so that
    the map epochs match the epochs the V2−V1 tables are evaluated at.
    """,
    source = @__FILE__,
    runs = [r.path for (r, _) in runs])
