"""
    BlockA

Parsers for the Block A community-benchmark runs, as they
are staged by `experiments/run_block_A.sh` under `runs/A/<id>/`.

The benchmark drivers in the model tree (`model/tests/test_benchmark_*.f90`)
write **only to stdout** — they read their vendored reference tables and print
fixed-width comparison tables plus a `PASS:`/`FAIL:` verdict. `runme` captures
that as `out.out`, together with the DKRZ epilogue that reports the job id,
elapsed time and peak RSS. So everything here is text parsing of `out.out`;
there is no NetCDF to read.

Alongside the text log, each driver now writes its comparison arrays — the
spatial profiles and time series behind every error number — as `#`-commented
column tables in the run directory (the opt-in `VILMA_BENCH_DUMP` path);
`dumptable` below reads one of those tables.

SIDE OF THE PIPELINE. This module reads `runs/` and is therefore part of the
REDUCTION half: it is included by `analysis/make_reduced.jl` and by nothing
else. The figures read `data/published/` through `analysis/Published.jl`, which
serves the same interface from the reduced files. If you are writing a figure
and reach for this module, you want Published instead.
"""
module BlockA

export Run, load_runs, verdict, passed, numeric_table, floats_in, slurm_meta, dumptable, hasdump, linear_nodes

const FLOAT_RE = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eEdD][-+]?\d+)?"

"""
Peak-RSS values at or below this are SLURM's unsampled floor, not measurements.

A step shorter than `JobAcctGatherFrequency` (30 s on levante) is never sampled,
and the epilogue then prints a constant -- 6100K -- for every such run. Several
block A benchmarks and the small-ℓmax block P runs finish inside that window, so
this has to be filtered or the figures report a 6 MB peak for a model that has
just read a 19 MB forcing file.
"""
const RSS_FLOOR_KB = 8192

"""
    Run

One benchmark run directory: its id (`"A1.love"`), the raw `out.out` lines, the
`PASS:`/`FAIL:` verdict line, and the SLURM accounting scraped from the DKRZ
epilogue.
"""
struct Run
    id::String
    dir::String
    lines::Vector{String}
    verdict::Union{String,Nothing}
    slurm::Dict{String,Any}
end

passed(r::Run) = r.verdict !== nothing && startswith(r.verdict, "PASS")

"""
    floats_in(line) -> Vector{Float64}

Every number on a line, regardless of the labels around it. Fortran `D`
exponents are normalised to `e`. Used for the label-and-unit-interleaved lines
(`"   homogeneous      -0.7306 m     -0.75 m     2.6 %"`).
"""
function floats_in(line::AbstractString)
    [parse(Float64, replace(m.match, r"[dD]" => "e")) for m in eachmatch(FLOAT_RE, line)]
end

"""
    numeric_table(lines, header; ncol, skiplabel=false) -> Matrix{Float64}

The block of numeric rows that follows the first line matching `header`.

Collection stops at the first line that is not a row of exactly `ncol` numbers,
so the surrounding prose terminates the table naturally. With
`skiplabel=true` the first whitespace-separated token is a string label
(`fig10`) and is dropped before parsing, so `ncol` counts the numbers only.

Returns a `nrow × ncol` matrix; empty if the header is not found.
"""
function numeric_table(lines::Vector{String}, header::Regex; ncol::Int, skiplabel::Bool=false)
    i = findfirst(l -> occursin(header, l), lines)
    i === nothing && return Matrix{Float64}(undef, 0, ncol)

    rows = Vector{Vector{Float64}}()
    for line in lines[(i+1):end]
        toks = split(strip(line))
        skiplabel && !isempty(toks) && (toks = toks[2:end])
        length(toks) == ncol || break
        vals = tryparse.(Float64, replace.(toks, r"[dD]" => "e"))
        any(isnothing, vals) && break
        push!(rows, Vector{Float64}(vals))
    end

    isempty(rows) && return Matrix{Float64}(undef, 0, ncol)
    permutedims(reduce(hcat, rows))
end

"""
    slurm_meta(lines) -> Dict

Job id, elapsed wall time (seconds), node and peak RSS (bytes) from the DKRZ
job epilogue appended to `out.out`. Missing keys simply do not appear — a run
that is still going has no epilogue yet.
"""
function slurm_meta(lines::Vector{String})
    meta = Dict{String,Any}()
    for l in lines
        if (m = match(r"^\*\s*JobID\s*:\s*(\d+)", l)) !== nothing
            meta["jobid"] = parse(Int, m[1])
        elseif (m = match(r"^\*\s*Elapsed time\s*:\s*(\d+):(\d\d):(\d\d)", l)) !== nothing
            meta["elapsed_s"] = 3600*parse(Int, m[1]) + 60*parse(Int, m[2]) + parse(Int, m[3])
        elseif (m = match(r"^\*\s*Nodelist\s*:\s*(\S+)", l)) !== nothing
            meta["node"] = m[1]
        elseif (m = match(r"^\*\s*Timelimit=(\d+):(\d\d):(\d\d)", l)) !== nothing ||
               (m = match(r"Timelimit=(\d+):(\d\d):(\d\d)", l)) !== nothing
            meta["wall_s"] = 3600*parse(Int, m[1]) + 60*parse(Int, m[2]) + parse(Int, m[3])
        elseif (m = match(r"^\*\s*0\s*\|\s*\S+\s+[\d.]+\s+(\d+)K", l)) !== nothing
            rss = parse(Int, m[1])
            # SLURM samples RSS on a timer (JobAcctGatherFrequency, 30 s here),
            # so a step that finishes inside one interval is never sampled and
            # the epilogue reports a fixed floor -- 6100K on levante, the same
            # value for every short run regardless of what it did. That is a
            # sentinel, not a measurement: the model cannot peak at 6 MB when it
            # has read a 19 MB forcing file. Recorded as missing rather than
            # plotted. See RSS_FLOOR_KB.
            rss > RSS_FLOOR_KB && (meta["maxrss_kb"] = rss)
        end
    end
    meta
end

"""
    load_runs(root) -> Vector{Run}

Every `<root>/<id>/out.out`, sorted by id. Directories without `out.out` (job
still queued) are skipped.
"""
function load_runs(root::AbstractString)
    runs = Run[]
    isdir(root) || return runs
    for id in sort(readdir(root))
        out = joinpath(root, id, "out.out")
        isfile(out) || continue
        lines = readlines(out)
        v = nothing
        for l in lines
            if (m = match(r"^\s*(PASS|FAIL):\s*(.*)$", l)) !== nothing
                v = string(m[1], ": ", strip(m[2]))
                break
            end
        end
        push!(runs, Run(id, joinpath(root, id), lines, v, slurm_meta(lines)))
    end
    runs
end

verdict(r::Run) = something(r.verdict, "(no verdict — job pending or crashed)")

"""
    dumptable(r, name) -> Matrix{Float64}

One of the driver's dumped column tables, e.g. `dumptable(run, "love_elastic.txt")`.
The leading `#` line names the columns and is skipped. Returns a 0×0 matrix if
the file is absent, so a caller can degrade gracefully when a run predates the
dump or the driver did not reach that point.
"""
function dumptable(r::Run, name::AbstractString)
    path = joinpath(r.dir, name)
    isfile(path) || return Matrix{Float64}(undef, 0, 0)
    rows = Vector{Vector{Float64}}()
    for line in eachline(path)
        startswith(lstrip(line), "#") && continue
        toks = split(strip(line))
        isempty(toks) && continue
        vals = tryparse.(Float64, toks)
        any(isnothing, vals) && continue
        push!(rows, Vector{Float64}(vals))
    end
    isempty(rows) && return Matrix{Float64}(undef, 0, 0)
    permutedims(reduce(hcat, rows))
end

"Column names from a dumped table's `#` header line."
function dumpcols(r::Run, name::AbstractString)
    path = joinpath(r.dir, name)
    isfile(path) || return String[]
    for line in eachline(path)
        startswith(lstrip(line), "#") && return split(strip(lstrip(lstrip(line), '#')))
    end
    String[]
end

hasdump(r::Run, name::AbstractString) = isfile(joinpath(r.dir, name))

"""
    linear_nodes(y; rtol=1e-9) -> Vector{Int}

Indices at which a piecewise-linear column carries independent data: the two
endpoints plus every interior point where the slope changes.

The Spada disc reference (`data/benchmarks/disc_spada2011/`) is distributed on a
0.1° θ grid but is piecewise linear on a ~1° one — all six time columns consist
of ~40 linear segments whose breaks fall at the same colatitudes, which no
smooth physical field can do (see that directory's PROVENANCE.md). Only the
breakpoints are data; the points between are straight-line fill. Comparing the
model against the full 0.1° grid therefore charges it for the reference's own
interpolation error, so the honest residual is evaluated at these indices.

The Martinec/SBK sea-level references are *not* like this — every point there
changes slope, i.e. every point is data — so this applies to A2 alone.

The threshold is **relative** to the largest curvature in the column. Within a
linear segment the second difference is exactly zero, and at a node it is O(the
curve's curvature), so the two are separated by many orders of magnitude and the
split is insensitive to `rtol`. An absolute threshold is not safe: it silently
turns into "every point is a node" the moment the values carry round-off, which
is exactly what a lower-precision dump caused here.
"""
function linear_nodes(y::AbstractVector; rtol = 1e-9)
    n = length(y)
    n < 3 && return collect(1:n)
    curv = abs.(diff(diff(y)))
    m = maximum(curv)
    m == 0 && return [1, n]                     # perfectly straight: no interior data
    brk = findall(>(rtol * m), curv) .+ 1       # curvature at index i sits at y[i+1]
    sort!(unique(vcat(1, brk, n)))
end

# --- per-benchmark extractors ------------------------------------------------
# Each returns a NamedTuple of whatever that driver prints. Fields are empty
# arrays when the driver did not get that far.

"A1 — `test_benchmark_love.x`: fluid-limit and elastic loading Love numbers."
function love(r::Run)
    fluid   = numeric_table(r.lines, r"h_ours\s+h_ref\s+l_ours"; ncol=7)
    elastic = numeric_table(r.lines, r"h_ours\s+h_ref\s+dh%";    ncol=7)
    (; fluid, elastic)
end

"A2 — `test_benchmark_disc.x`: profile error summary + VE centre uplift u(0,t)."
function disc(r::Run)
    centre = numeric_table(r.lines, r"t\[kyr\]\s+mine\[m\]\s+Spada\[m\]"; ncol=4)
    uline = findfirst(l -> occursin(r"^\s+u: peak", l), r.lines)
    nline = findfirst(l -> occursin(r"^\s+N: peak", l), r.lines)
    u = uline === nothing ? Float64[] : floats_in(r.lines[uline])   # peak, max|du|, %, theta
    n = nline === nothing ? Float64[] : floats_in(r.lines[nline])   # peak, max|dN|, %
    (; centre, u_summary=u, n_summary=n)
end

"A3 case A — `test_benchmark_martinec.x`: peak and max relative error per field."
function martinec_caseA(r::Run)
    fields = String[]; peaks = Float64[]; errs = Float64[]
    for l in r.lines
        m = match(r"^\s+(uplift|horizontal|geoid)\s+peak\s", l)
        m === nothing && continue
        v = floats_in(l)
        length(v) < 2 && continue
        push!(fields, m[1]); push!(peaks, v[1]); push!(errs, v[2])
    end
    (; fields, peaks, errs)
end

"""
A3 C2/D3/E2/F1 — `test_benchmark_sle.x`: the peak-normalized max error table
(one row per Martinec figure, columns u / vth / vph / N / ss / sle), plus the
eustatic and mass-conservation diagnostics.
"""
function sle(r::Run)
    tab  = numeric_table(r.lines, r"fig\s+rows\s+u\[%\]"; ncol=7, skiplabel=true)
    figs = String[]
    i = findfirst(l -> occursin(r"fig\s+rows\s+u\[%\]", l), r.lines)
    if i !== nothing
        for l in r.lines[(i+1):end]
            m = match(r"^\s+(fig\d+)\s", l)
            m === nothing && break
            push!(figs, m[1])
        end
    end
    eus  = findfirst(l -> occursin(r"^\s*eustatic:", l), r.lines)
    ofrac = findfirst(l -> occursin(r"^\s*ocean frac:", l), r.lines)
    tally = findfirst(l -> occursin(r"SLE iteration tally", l), r.lines)
    (; figs, table=tab,
       eustatic = eus === nothing ? Float64[] : floats_in(r.lines[eus]),
       ocean    = ofrac === nothing ? Float64[] : floats_in(r.lines[ofrac]),
       iters    = tally === nothing ? Float64[] : floats_in(r.lines[tally]))
end

"""
A4 — `test_rotation.x`: polar motion |m|(t) against Spada Table 14, for both
the cap and the disc load, plus the degree-2 tidal Love numbers.
"""
function rotation(r::Run)
    bodies = Dict{String,Matrix{Float64}}()
    G = Dict{String,Float64}()
    for (i, l) in enumerate(r.lines)
        m = match(r"^\s+\(3\)\s+(\w+): polar motion", l)
        if m !== nothing
            bodies[m[1]] = numeric_table(r.lines[(i+1):end], r"t\[kyr\]\s+\|m\| model"; ncol=4)
        end
        m = match(r"^\s+\(1\)\s+(\w+): rigid", l)
        if m !== nothing
            # the rel.err line follows two lines later
            j = findfirst(k -> occursin(r"rel\.err\s*=", r.lines[k]), (i+1):min(i+5, length(r.lines)))
            j !== nothing && (G[m[1]] = floats_in(r.lines[i+j])[1])
        end
    end
    kte  = findfirst(l -> occursin(r"elastic\s+k\^T_e", l), r.lines)
    ks   = findfirst(l -> occursin(r"secular\s+k_s", l), r.lines)
    (; bodies, G,
       kTe = kte === nothing ? NaN : floats_in(r.lines[kte])[end],
       k_s = ks  === nothing ? NaN : floats_in(r.lines[ks])[end])
end

"""
A5 — `test_benchmark_lvz.x`: centre uplift U(t) for the LVZ case, and the
200-yr centre values for the homogeneous and LVZ Earths against Weerdesteijn.
"""
function lvz(r::Run)
    t = Float64[]; u = Float64[]
    for l in r.lines
        m = match(r"^\s+t=\s*([\d.]+)\s*yr\s+U_center=\s*(\S+)", l)
        m === nothing && continue
        push!(t, parse(Float64, m[1])); push!(u, parse(Float64, m[2]))
    end
    cases = Dict{String,NTuple{3,Float64}}()   # model, reference, rel.err [%]
    for l in r.lines
        m = match(r"^\s+(homogeneous|with LVZ)\s", l)
        m === nothing && continue
        v = floats_in(l)
        length(v) >= 3 && (cases[m[1]] = (v[1], v[2], v[3]))
    end
    amp = findfirst(l -> occursin(r"LVZ amplification", l), r.lines)
    (; t, u, cases, amplification = amp === nothing ? Float64[] : floats_in(r.lines[amp]))
end

end # module
