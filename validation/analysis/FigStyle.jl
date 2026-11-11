"""
    FigStyle

What a figure of the validation suite looks like, what it is called and where it
goes.

  * A figure carries axes, units, a legend and panel labels (a), (b), …; the
    description is its caption, written beside it as `figures/<name>.md`.
  * The figures the documentation shows are listed in `DOCS_FIGURES`, slug =>
    file name, and are written to `docs/benchmarks/figures/` (tracked). Every
    other figure is a diagnostic and goes to `validation/figures/` (untracked),
    as do all captions.

    using .FigStyle
    fig = Figure(size = (1100, 620))
    panel!(fig, (1, 1), "a")
    savefig(fig, "love-numbers"; caption = "...", source = @__FILE__)

Scripts also `record!` their headline numbers into `tables/results.toml`.
"""
module FigStyle

using CairoMakie
using Printf
using Statistics
using TOML

export DOCS_FIGURES, panel!, panels!, savefig, record!, results_path, tables_dir,
       figures_dir, docs_figures_dir, C_1D, C_3D, C_TOR, C_V2, C_V1, C_GREY,
       CMAP_DIV, CMAP_SEQ, symlim, TICK_FS,
       SCALE_RUNS, SCALE_EPOCHS, SCALE_DIFF_PAIRS, deglac_scales

const REPO  = dirname(@__DIR__)        # validation/
const MODEL = dirname(REPO)            # the repository root

# The documentation set: slug passed to `savefig` => file name (no extension)
# under docs/benchmarks/figures/.
const DOCS_FIGURES = Dict(
    # block A — community benchmarks
    "love-numbers"   => "love-numbers",
    "disc-benchmark" => "disc",
    "sle-benchmark"  => "sle",
    "sle-profiles-C2" => "sle-C2",
    "sle-profiles-D3" => "sle-D3",
    "sle-profiles-E2" => "sle-E2",
    "sle-profiles-F1" => "sle-F1",
    "polar-motion"   => "rotation",
    "lvz-benchmark"  => "lvz",
    # block B — lateral-viscosity benchmark
    "bench-visc3d-sections-heav"  => "visc3d-sections-heav",
    "bench-visc3d-sections-ramp"  => "visc3d-sections-ramp",
    "bench-visc3d-geocentre"      => "visc3d-geocentre",
    "bench-visc3d-A-normal-mode"  => "visc3d-A-normal-mode",
    "bench-visc3d-lmax"           => "visc3d-lmax",
    # last deglaciation — VILMA2 and VILMA2 against VILMA1
    "deglaciation-maps" => "deglac-maps",
    "v2v1-maps-3d"      => "deglac-v2v1-maps-3d",
    "v2v1-maps-1db"     => "deglac-v2v1-maps-1d",
    "v2v1-summary"      => "deglac-v2v1-summary",
)

figures_dir()      = joinpath(REPO, "figures")
docs_figures_dir() = joinpath(MODEL, "docs", "benchmarks", "figures")
tables_dir()       = mkpath(joinpath(REPO, "tables"))
results_path()     = joinpath(tables_dir(), "results.toml")

# ---------------------------------------------------------------------------
# Palette and theme
# ---------------------------------------------------------------------------
# Two axes of meaning, kept apart so a figure never has to explain its colours:
# Earth structure (1-D vs 3-D) and solver (VILMA2 vs VILMA1).
const C_1D   = "#1f5c99"
const C_3D   = "#c8553d"
const C_TOR  = "#6e1f12"     # 3-D with the toroidal coupling: a darker 3-D red
const C_V2    = "#c8553d"     # VILMA2 — the same red as 3-D: V2 is the 3-D code
const C_V1    = "#2f6f4e"
const C_GREY = "#8a8f98"

const CMAP_DIV = :balance                                    # signed fields
const CMAP_SEQ = cgrad([:grey25, :grey55, :grey85, :white])  # ice thickness overlay

# Point sizes are set for a slide first and the page second. TICK_FS is shared
# with the map scripts, which override the axis tick sizes.
const TICK_FS = 12

const THEME = Theme(
    fontsize = 14,
    figure_padding = 10,
    Axis = (
        titlesize = 14, titlefont = :regular, titlealign = :left,
        xlabelsize = 14, ylabelsize = 14,
        xticklabelsize = TICK_FS, yticklabelsize = TICK_FS,
        xgridvisible = false, ygridvisible = false,
        spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8,
    ),
    Legend = (framevisible = false, labelsize = 12, patchsize = (20, 12)),
    Colorbar = (labelsize = 13, ticklabelsize = 12, size = 11),
    Lines = (linewidth = 2.0,),
    Scatter = (markersize = 8, strokewidth = 0),
)

function __init__()
    set_theme!(THEME)
    CairoMakie.activate!(type = "png", px_per_unit = 2)
end

# ---------------------------------------------------------------------------
# Panel labels
# ---------------------------------------------------------------------------
"""
    panel!(fig, pos, letter)

Put `(a)`, `(b)`, … in the top-left gutter of the layout cell at `pos`, outside
the axis so it never lands on data. `pos` is whatever indexes the layout:
`(1, 2)`, `(1, 1:3)`, …
"""
function panel!(fig, pos, letter::AbstractString)
    Label(fig[pos..., TopLeft()], "($letter)";
          font = :bold, fontsize = 15, halign = :left,
          padding = (0, 0, 2, 0))
end
panel!(fig, pos, i::Integer) = panel!(fig, pos, string('a' + i - 1))

"""
    panels!(fig, positions)

Label a whole figure at once, in the order given: `panels!(fig, [(1,1), (1,2), (2,1)])`.
"""
panels!(fig, positions) = [panel!(fig, p, i) for (i, p) in enumerate(positions)]

# ---------------------------------------------------------------------------
# Colour limits
# ---------------------------------------------------------------------------
"""
    symlim(fields...; q = 0.995) -> (-m, m)

Symmetric colour limits at the `q` quantile of |value|, so that a handful of
extreme cells — which on a GIA field means the middle of Hudson Bay — do not
flatten the rest of the map to white.
"""
function symlim(fields...; q = 0.995)
    v = Float64[]
    for f in fields
        append!(v, filter(isfinite, vec(f)))
    end
    isempty(v) && return (-1.0, 1.0)
    m = quantile!(abs.(v), q)
    m = m > 0 ? m : maximum(abs.(v); init = 1.0)
    (-m, m)
end

# ---------------------------------------------------------------------------
# Colour scales shared BETWEEN figures
# ---------------------------------------------------------------------------
# A colour range computed per figure makes every difference map look equally
# large. The deglaciation maps therefore take their limits from one place,
# computed over all the runs and epochs the set draws, so that a panel is
# comparable with any other panel of the same kind in any other figure.

"The runs whose fields share a colour scale: the block D production matrix."
const SCALE_RUNS = ["D-V2-3D.l170.fe", "D-V2-1Db.l170.fe",
                    "D-V1-3D.l170.fe", "D-V1-1Db.l170.fe"]

"The epochs any block D map is drawn at. All must be in `Published.EPOCHS`."
const SCALE_EPOCHS = [-21000.0, -14000.0, -10000.0, -8000.0]

"""
Every model-to-model difference the figure set draws, as `(a, b)` run ids
meaning `a − b`. The lateral-viscosity signal and the two solver residuals are
deliberately in one list: they are the same KIND of quantity — one model's field
minus another's on the identical problem — and the whole point is that a reader
can weigh them against each other.
"""
const SCALE_DIFF_PAIRS = [
    ("D-V2-3D.l170.fe",  "D-V2-1Db.l170.fe"),   # lateral structure, solver held fixed
    ("D-V2-3D.l170.fe",  "D-V1-3D.l170.fe"),    # solver, 3-D structure held fixed
    ("D-V2-1Db.l170.fe", "D-V1-1Db.l170.fe"),   # solver, 1Db structure held fixed
]

"""
    deglac_scales(runs, rsl_anomaly, uplift_rate) -> (; rsl, upl, drsl, dupl)

The four colour ranges every block D map uses:

    rsl   any panel showing relative sea level vs present day      [m]
    upl   any panel showing present-day uplift rate                [mm/yr]
    drsl  any panel showing one model's rsl MINUS another's        [m]
    dupl  the same for uplift rate                                 [mm/yr]

`runs` is a `Dict` of id => run covering `SCALE_RUNS`; the two accessors are
passed in rather than imported so that this module stays independent of how the
data is loaded. Every caller must pass the SAME dict — which is why the three
map scripts each load the full production matrix even when they plot a subset of
it. Load only the runs you draw and you are back to a per-figure scale.

Signal and difference keep separate ranges on purpose. On the signal range a
residual of a few per cent renders as an empty white panel, which would show
agreement rather than demonstrate it; the split is what makes the difference
columns readable at all. What is NOT allowed is a different difference range per
figure.

A run that was not published simply drops out of the reduction, which changes the
scale — so a partial set is internally consistent but not comparable with a full
one. `savefig` records which runs a figure came from, and that is where to look
if two figures will not line up.
"""
function deglac_scales(runs::AbstractDict, rsl_anomaly, uplift_rate)
    have(id) = haskey(runs, id) && runs[id] !== nothing
    ids = filter(have, SCALE_RUNS)
    isempty(ids) && error("deglac_scales: none of $(SCALE_RUNS) is loaded")
    length(ids) == length(SCALE_RUNS) ||
        @warn "deglac_scales: only $(length(ids)) of $(length(SCALE_RUNS)) runs loaded — " *
              "the colour scale is not comparable with a figure made from the full set" ids

    sig_rsl = [rsl_anomaly(runs[i], t) for i in ids for t in SCALE_EPOCHS]
    sig_upl = [uplift_rate(runs[i]) for i in ids]

    pairs = [p for p in SCALE_DIFF_PAIRS if have(p[1]) && have(p[2])]
    dif_rsl = [rsl_anomaly(runs[a], t) .- rsl_anomaly(runs[b], t)
               for (a, b) in pairs for t in SCALE_EPOCHS]
    dif_upl = [uplift_rate(runs[a]) .- uplift_rate(runs[b]) for (a, b) in pairs]

    (rsl  = symlim(sig_rsl...),
     # q = 0.999 on the uplift rate: it is concentrated in a few hundred cells
     # over Hudson Bay and Fennoscandia, and the 0.995 quantile clips the peaks
     # that are the whole point of the panel.
     upl  = symlim(sig_upl...; q = 0.999),
     drsl = isempty(dif_rsl) ? (-1.0, 1.0) : symlim(dif_rsl...),
     dupl = isempty(dif_upl) ? (-1.0, 1.0) : symlim(dif_upl...; q = 0.999))
end

# ---------------------------------------------------------------------------
# Saving, with the caption
# ---------------------------------------------------------------------------
"""
    savefig(fig, slug; caption, title, notes, source, runs) -> path

Write the PNG (to `docs/benchmarks/figures/` under its `DOCS_FIGURES` name, else
to `validation/figures/<slug>.png`) and the caption `figures/<name>.md`.
`caption` describes the figure completely, naming every panel by its letter;
`notes` says what the figure is evidence for. `runs` are recorded for provenance.
"""
function savefig(fig, slug::AbstractString;
                 caption::AbstractString,
                 title::AbstractString = "",
                 notes::AbstractString = "",
                 source::AbstractString = "",
                 runs::AbstractVector = String[])
    name = get(DOCS_FIGURES, String(slug), String(slug))
    dir  = haskey(DOCS_FIGURES, String(slug)) ? docs_figures_dir() : figures_dir()
    mkpath(dir); mkpath(figures_dir())
    png = joinpath(dir, name * ".png")
    save(png, fig)

    open(joinpath(figures_dir(), name * ".md"), "w") do io
        println(io, "# ", isempty(title) ? name : title, "\n")
        println(io, strip(caption), "\n")
        if !isempty(strip(notes))
            println(io, "**Notes.** ", strip(notes), "\n")
        end
        src = isempty(source) ? "" : relpath(source, MODEL)
        prov = String[]
        isempty(src) || push!(prov, "`$src`")
        isempty(runs) || push!(prov, join(["`" * relpath(String(r), MODEL) * "`" for r in runs], ", "))
        isempty(prov) || println(io, "*Generated by ", join(prov, " from "), ".*")
    end
    @info "wrote $(relpath(png, MODEL))"
    png
end

# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------
"""
    record!(block, key, value)

Store one headline number (or string, or vector) under `[block]` in
`tables/results.toml`. Read-modify-write on every call, so the scripts can run
individually and in any order, and a rerun of one block updates only its keys.
"""
function record!(block::AbstractString, key::AbstractString, value)
    p = results_path()
    d = isfile(p) ? TOML.parsefile(p) : Dict{String,Any}()
    b = get!(d, block, Dict{String,Any}())
    b[key] = value
    mkpath(dirname(p))
    open(p, "w") do io
        println(io, "# Headline numbers, written by validation/analysis/*.jl via FigStyle.record!.")
        println(io, "# Generated file — do not edit by hand.\n")
        TOML.print(io, d; sorted = true)
    end
    value
end

"""
    record!(block, pairs...)

`record!("blockD", "bsl_3d" => 98.32, "bsl_1db" => 98.36)` — one file rewrite
for the lot.
"""
function record!(block::AbstractString, pairs::Pair...)
    p = results_path()
    d = isfile(p) ? TOML.parsefile(p) : Dict{String,Any}()
    b = get!(d, block, Dict{String,Any}())
    for (k, v) in pairs
        b[String(k)] = v
    end
    mkpath(dirname(p))
    open(p, "w") do io
        println(io, "# Headline numbers, written by validation/analysis/*.jl via FigStyle.record!.")
        println(io, "# Generated file — do not edit by hand.\n")
        TOML.print(io, d; sorted = true)
    end
    d
end

end # module
