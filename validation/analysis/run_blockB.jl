#!/usr/bin/env julia
#
# run_blockB.jl — block B, the lateral-viscosity benchmark (Klemann et al. setup
# note; model docs/benchmarks/visc3d.qmd, driver tests/bench_visc3d.f90).
#
#   julia --project=. analysis/run_blockB.jl
#
# Tests A-i (no structure), B-i (low-viscosity column under the load, at the
# pole) and C-i (column off the load centre, fully 3-D), each under a Heaviside
# and a 10-kyr ramp load, ℓmax 128, toroidal coupling on, to 100 kyr. The
# benchmark provides no reference solution, so this block shows VILMA2's answer
# in the protocol's own form and records the numbers a comparison will be read on.
#
# Two checks it can make on its own:
#   * A-i is laterally uniform, so its exact answer is the normal-mode solution
#     of M3-L70-V01 (analysis/NormalMode.jl, from the published Love table).
#   * A-i and B-i are axisymmetric and cheap, so they are also run at ℓmax 64
#     and 256: A against the normal-mode solution to degree 256, B against its
#     own ℓmax 256 run.
#
# Drawn from data/published/bench_visc3d/, the driver's protocol files as
# submitted. The layout follows the model's analysis/plot_bench_visc3d.jl.

using CairoMakie
using Printf

include(joinpath(@__DIR__, "Published.jl"))
include(joinpath(@__DIR__, "FigStyle.jl"))
include(joinpath(@__DIR__, "NormalMode.jl"))
using .Published
using .FigStyle
using .NormalMode

const REPO = dirname(@__DIR__)

cases = load_blockB()
isempty(cases) && error("no block B case is published — run ./run.sh reduce")
# The convergence truncations; an ℓmax that was not published is simply absent.
conv = Dict(L => load_blockB(L) for L in BENCH_CONV_LMAX if isdir(Published.bench_case_dir(L)))
conv[BENCH_LMAX] = cases

love = read_love_table(bench_love_path())
σcap = cap_legendre(love.nmax)       # the A-i load, degrees 0..256

# Protocol columns (Published.BenchV3D).
const COL_T = 4
const FIELDS = [(5, L"u_r\ \mathrm{[m]}"), (6, L"u_\vartheta\ \mathrm{[m]}"),
                (7, L"u_\varphi\ \mathrm{[m]}"), (8, L"\delta\phi\ \mathrm{[m^2\,s^{-2}]}")]

# Distance of the load centre from the structure centre along the section [deg]:
# A and B put both at the pole; C puts the structure at colatitude 35° and the
# load at 30° on the same meridian (bench_visc3d.f90).
const LOAD_DIST = Dict("A" => 0.0, "B" => 0.0, "C" => 5.0)

const FORCING_LABEL = Dict("heav" => "Heaviside load", "ramp" => "load ramped over 0–10 kyr")
const TRANGE = (-1.0, 2.0)                  # log10 of the epoch range [kyr]

lmax_of(cs) = only(unique(c.lmax for c in values(cs)))

"Rows of one epoch."
epoch_rows(c, tk) = c.disp[:, COL_T] .== tk

"The section value of column `col` nearest the load centre, at epoch `tk`."
function at_load(c, col, tk)
    d = c.disp[epoch_rows(c, tk), :]
    d[argmin(abs.(d[:, 3] .- LOAD_DIST[c.test])), col]
end

# =============================================================================
# Cross-sections, one figure per forcing
# =============================================================================
function fig_sections(forcing)
    have = [t for t in BENCH_TESTS if haskey(cases, (t, forcing))]
    isempty(have) && (@warn "no $forcing case published"; return nothing)
    fig = Figure(size = (430 * length(have) + 160, 1050))
    pos = Tuple{Int,Int}[]
    for (j, t) in enumerate(have), (i, (col, lab)) in enumerate(FIELDS)
        ax = Axis(fig[i, j]; ylabel = j == 1 ? lab : "",
                  xlabel = i == length(FIELDS) ? "distance from structure centre [deg]" : "",
                  title = i == 1 ? "$(t)-i" : "",
                  xticks = 0:30:180)
        xlims!(ax, 0, 180)
        push!(pos, (i, j))
        haskey(cases, (t, forcing)) || continue
        c = cases[(t, forcing)]
        for tk in BENCH_EPOCHS_KYR
            sel = epoch_rows(c, tk)
            lines!(ax, c.disp[sel, 3], c.disp[sel, col]; color = log10(tk),
                   colorrange = TRANGE, colormap = :viridis, linewidth = 1.5)
        end
        # An identically zero field gets Makie's default ±1 range, whose top tick
        # label sits on the axis top and collides with the panel letter.
        maximum(abs.(c.disp[:, col])) < 1e-6 && ylims!(ax, -1.25, 1.25)
    end
    Colorbar(fig[:, length(have) + 1]; colormap = :viridis, limits = TRANGE,
             label = "time [kyr]", ticks = (-1:2, ["0.1", "1", "10", "100"]))
    colgap!(fig.layout, 44)
    # Label down the columns, then across: (a)-(d) are test A's four fields.
    panels!(fig, sort(pos, by = p -> (p[2], p[1])))
    savefig(fig, "bench-visc3d-sections-$forcing";
        title = "Lateral-viscosity benchmark: cross-sections, $(FORCING_LABEL[forcing])",
        caption = """
        VILMA2's solution of the lateral-viscosity benchmark of Klemann et al., tests
        $(join(["$(t)-i" for t in have], ", ", " and ")) (columns), $(lowercase(FORCING_LABEL[forcing])), at ℓmax
        $(lmax_of(cases)). Rows: radial displacement u_r, colatitudinal and
        longitudinal displacement u_ϑ and u_φ, and the perturbation of the
        gravitational potential δφ at the surface, along the section from the centre
        of the low-viscosity structure through the load centre to the structure's
        antipode. Colour is time, at the protocol's thirteen epochs from 0.1 to
        100 kyr. (a–d) Test A-i, M3-L70-V01 with no structure; (e–h) test B-i, a
        column of 12° half-width from 60 to 220 km depth with viscosity reduced to
        10¹⁸ Pa s at its centre, under the load at the pole$("C" in have ? "; (i–l) test C-i, the same
        column centred 5° from the load centre, which makes the problem fully 3-D" : ""). $("C" in have ? "" :
        "Test C-i (the column centred 5° from the load centre, fully 3-D) was not run in this pass.")
        """,
        notes = """
        There is no reference solution: the benchmark collects submissions from
        several codes and compares them, so this figure is VILMA2's entry in the
        protocol's own form, from the files submitted
        (`data/published/bench_visc3d/`). u_φ is zero on the section in every test:
        the section runs along the meridian through both centres, a plane of mirror
        symmetry of the problem, and what C-i shows (~10⁻⁶ m) is round-off. The jump
        in C-i's u_ϑ at 35° is where the section crosses the pole: the colatitude
        direction reverses there, so the component changes sign by definition. The
        ramped load reaches its full thickness at 10 kyr, which is why the 11, 12
        and 15 kyr epochs are in the protocol. The toroidal coupling is on
        (`l_toroidal`), as the lateral structure drives a toroidal flow.
        """,
        source = @__FILE__,
        runs = [joinpath(REPO, "data", "published", "bench_visc3d")])
end

# =============================================================================
# Geocentre motion
# =============================================================================
function fig_gcm()
    fig = Figure(size = (1250, 620))
    colors = Dict("A" => C_GREY, "B" => C_1D, "C" => C_3D)
    pos = Tuple{Int,Int}[]
    for (i, f) in enumerate(BENCH_FORCINGS), (j, comp) in enumerate(["u_x", "u_y", "u_z"])
        ax = Axis(fig[i, j]; xscale = log10, xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"]),
                  xlabel = i == length(BENCH_FORCINGS) ? "time [kyr]" : "",
                  ylabel = j == 1 ? "$(FORCING_LABEL[f])\n[m]" : "",
                  title = i == 1 ? comp : "")
        push!(pos, (i, j))
        for t in BENCH_TESTS
            haskey(cases, (t, f)) || continue
            g = cases[(t, f)].gcm
            scatterlines!(ax, g[:, 1], g[:, j + 1]; color = colors[t], label = "$(t)-i",
                          markersize = 6)
        end
        # zero components (A-i and B-i are axisymmetric): keep the top tick off the panel letter
        all(maximum(abs.(cases[(t, f)].gcm[:, j + 1])) < 1e-6
            for t in BENCH_TESTS if haskey(cases, (t, f))) && ylims!(ax, -1.25, 1.25)
        (i, j) == (1, 3) && axislegend(ax; position = :lb)
    end
    colgap!(fig.layout, 44)
    panels!(fig, pos)
    savefig(fig, "bench-visc3d-geocentre";
        title = "Lateral-viscosity benchmark: geocentre motion",
        caption = """
        Geocentre motion u_CF − u_CM in the benchmark of Klemann et al., Cartesian
        components (columns) against time, for tests $(join(["$(t)-i" for t in BENCH_TESTS if any(haskey(cases, (t, f)) for f in BENCH_FORCINGS)], ", ", " and ")), under (a–c) the
        Heaviside load and (d–f) the load ramped over 0–10 kyr, at ℓmax
        $(lmax_of(cases)). The model runs in the centre-of-mass frame, so u_CF − u_CM
        is the surface mean of the displacement vector.
        """,
        notes = """
        A-i and B-i put the load and the structure at the pole, so only u_z is
        non-zero there. In C-i both sit on the 25°E meridian at colatitudes 30° (load)
        and 35° (structure), which gives x and y components as well.
        """,
        source = @__FILE__,
        runs = [joinpath(REPO, "data", "published", "bench_visc3d")])
end

# =============================================================================
# A-i against the normal-mode solution
# =============================================================================
"Normal-mode u_r, u_ϑ, δφ on case `c`'s section at epoch `tk`, summed to degree `nmax`."
function nm_section(c, tk; nmax = c.lmax)
    d = c.disp[epoch_rows(c, tk), 3]
    section_ref(love, σcap, d, tk, c.forcing; nmax)
end

"max |VILMA2 − normal mode| of u_r, u_ϑ, δφ on the section at `tk` (same truncation)."
function nm_error(c, tk)
    ur, ut, dp = nm_section(c, tk)
    d = c.disp[epoch_rows(c, tk), :]
    (ur = maximum(abs.(d[:, 5] .- ur)), ut = maximum(abs.(d[:, 6] .- ut)),
     dp = maximum(abs.(d[:, 8] .- dp)), ur_nm0 = ur[1], ur0 = d[1, 5],
     ur_scale = maximum(abs.(ur)))
end

function fig_normal_mode()
    have = [f for f in BENCH_FORCINGS if haskey(cases, ("A", f))]
    isempty(have) && (@warn "A-i not published"; return nothing)
    fig = Figure(size = (1150, 1150))
    rows = [(5, L"\Delta u_r\ \mathrm{[m]}", 1), (6, L"\Delta u_\vartheta\ \mathrm{[m]}", 2),
            (8, L"\Delta\delta\phi\ \mathrm{[m^2\,s^{-2}]}", 3)]
    pos = Tuple{Int,Int}[]
    for (j, f) in enumerate(BENCH_FORCINGS)
        haskey(cases, ("A", f)) || continue
        c = cases[("A", f)]
        for (i, (col, lab, k)) in enumerate(rows)
            ax = Axis(fig[i, j]; ylabel = j == 1 ? lab : "", xticks = 0:30:180,
                      title = i == 1 ? uppercasefirst(FORCING_LABEL[f]) : "")
            xlims!(ax, 0, 180); push!(pos, (i, j))
            for tk in BENCH_EPOCHS_KYR
                ref = nm_section(c, tk)[k]
                d = c.disp[epoch_rows(c, tk), :]
                lines!(ax, d[:, 3], d[:, col] .- ref; color = log10(tk), colorrange = TRANGE,
                       colormap = :viridis, linewidth = 1.5)
            end
        end
        ax = Axis(fig[4, j]; xscale = log10, xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"]),
                  xlabel = "time [kyr]", ylabel = j == 1 ? L"u_z\ \mathrm{geocentre\ [m]}" : "")
        push!(pos, (4, j))
        tt = 10 .^ range(-1, 2, length = 300)
        lines!(ax, tt, [gcm_ref(love, σcap, t, f) for t in tt]; color = C_GREY, label = "normal modes")
        scatter!(ax, c.gcm[:, 1], c.gcm[:, 4]; color = C_3D, label = "VILMA2")
        j == 1 && axislegend(ax; position = :lb)
    end
    Label(fig[3, 1:2, Bottom()], "distance from the load centre [deg]"; padding = (0, 0, 0, 30))
    Colorbar(fig[1:3, 3]; colormap = :viridis, limits = TRANGE, label = "time [kyr]",
             ticks = (-1:2, ["0.1", "1", "10", "100"]))
    rowgap!(fig.layout, 3, 45)
    panels!(fig, sort(pos, by = p -> (p[2], p[1])))
    e = [nm_error(cases[("A", f)], tk) for f in have for tk in BENCH_EPOCHS_KYR]
    savefig(fig, "bench-visc3d-A-normal-mode";
        title = "Lateral-viscosity benchmark: test A-i against the normal-mode solution",
        caption = """
        Test A-i of the lateral-viscosity benchmark, where the Earth is laterally uniform
        M3-L70-V01 and the exact answer is known: VILMA2 minus the normal-mode solution,
        summed to the same truncation ℓmax $(lmax_of(cases)), along the section from the
        load centre, at the protocol's thirteen epochs (colour). (a–c) Heaviside load:
        radial and colatitudinal displacement and the surface potential perturbation δφ;
        (d) the geocentre motion u_CF − u_CM along the load axis, VILMA2 (points) and
        normal modes (line). (e–h) The same for the load ramped over 0–10 kyr. Over all
        epochs and both load histories the largest differences are
        $(@sprintf("%.3f", maximum(x.ur for x in e))) m in u_r (against a signal of
        $(@sprintf("%.0f", maximum(x.ur_scale for x in e))) m),
        $(@sprintf("%.3f", maximum(x.ut for x in e))) m in u_ϑ and
        $(@sprintf("%.3f", maximum(x.dp for x in e))) m² s⁻² in δφ.
        """,
        notes = """
        The normal-mode Love numbers are giapy's table for M3-L70-V01 (Spada et al. 2011,
        incompressible, nine modes per degree, degrees 1–256), vendored in the model repo
        and published beside the runs. The benchmark's Earth differs from it only by two
        added interfaces (60, 220 km) and a 60–70 km layer of 10³⁰ Pa s instead of
        elastic, neither of which matters over 100 kyr. Degree 1 is in the
        centre-of-mass frame in both. The ramp is the analytic convolution of the
        Heaviside response with the linear load history; the load's own potential in δφ
        follows the load fraction directly. `analysis/NormalMode.jl` has the formulas.
        """,
        source = @__FILE__,
        runs = [joinpath(REPO, "data", "published", "bench_visc3d")])
end

# =============================================================================
# Truncation: ℓmax 64 / 128 / 256 for A-i and B-i
# =============================================================================
"rms over the section of u_r(a) − u_r(b) at `tk` (b is a field on the same points)."
rms_diff(a, b) = sqrt(sum(abs2, a .- b) / length(a))

function conv_series()
    out = Dict{Tuple{String,String,Int},Vector{Float64}}()
    Ls = sort(collect(keys(conv)))
    for f in BENCH_FORCINGS, L in Ls
        cA = get(conv[L], ("A", f), nothing)
        if cA !== nothing
            out[("A", f, L)] = [rms_diff(cA.disp[epoch_rows(cA, tk), 5],
                                         nm_section(cA, tk; nmax = love.nmax)[1])
                                for tk in BENCH_EPOCHS_KYR]
        end
        cB, cR = get(conv[L], ("B", f), nothing), get(get(conv, 256, Dict()), ("B", f), nothing)
        if cB !== nothing && cR !== nothing && L != 256
            out[("B", f, L)] = [rms_diff(cB.disp[epoch_rows(cB, tk), 5], cR.disp[epoch_rows(cR, tk), 5])
                                for tk in BENCH_EPOCHS_KYR]
        end
    end
    out
end

const SUPERSCRIPT = Dict(zip("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
"Tick label of 10^d: `1`, `10` for d ≥ 0, else `10⁻³`."
decade_label(d) = d >= 0 ? string(10^d) : "10" * map(c -> SUPERSCRIPT[c], string(d))

function fig_lmax(cs)
    isempty(cs) && (@warn "no convergence run published"; return nothing)
    fig = Figure(size = (1100, 700))
    col = Dict(64 => C_1D, 128 => C_3D, 256 => C_TOR)
    pos = Tuple{Int,Int}[]
    for (i, t) in enumerate(["A", "B"]), (j, f) in enumerate(BENCH_FORCINGS)
        ref = t == "A" ? "normal modes to degree $(love.nmax)" : "B-i at ℓmax 256"
        # Decade ticks only: Makie labels a sub-decade range as 10^0.5, which is not
        # a number a reader converts in their head (as in run_blockE.jl).
        dec = -6:1
        ax = Axis(fig[i, j]; xscale = log10, yscale = log10,
                  xticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"]),
                  yticks = (10.0 .^ dec, decade_label.(dec)),
                  xlabel = i == 2 ? "time [kyr]" : "",
                  ylabel = j == 1 ? "$(t)-i: rms Δu_r [m]\nvs $ref" : "",
                  title = i == 1 ? uppercasefirst(FORCING_LABEL[f]) : "")
        push!(pos, (i, j))
        for L in (64, 128, 256)
            haskey(cs, (t, f, L)) || continue
            v = cs[(t, f, L)]
            keep = v .> 0
            scatterlines!(ax, BENCH_EPOCHS_KYR[keep], v[keep]; color = col[L], label = "ℓmax $L",
                          markersize = 6)
        end
        (i, j) == (1, 1) && axislegend(ax; position = :rb)
    end
    panels!(fig, pos)
    savefig(fig, "bench-visc3d-lmax";
        title = "Lateral-viscosity benchmark: truncation",
        caption = """
        Truncation of the axisymmetric benchmark tests at ℓmax 64, 128 and 256: root-mean-
        square difference in u_r along the 0–180° section. (a, b) A-i against the
        normal-mode solution summed to degree $(love.nmax), which measures truncation and
        solver error together; (c, d) B-i, which has no exact solution, against its own
        ℓmax 256 run. Left, Heaviside load; right, load ramped over 0–10 kyr.
        """,
        notes = """
        The rms is taken over the section's points, spaced uniformly in distance, not
        weighted by area. For A-i at ℓmax 256 the difference is the solver error alone,
        since the table stops at degree 256. The model docs report Gibbs ringing at ℓmax
        64 at the load edge, which is where the ℓmax 64 difference sits.
        """,
        source = @__FILE__,
        runs = [joinpath(REPO, "data", "published", "bench_visc3d")])
end

# =============================================================================
# Drive, and the table
# =============================================================================
for f in BENCH_FORCINGS
    fig_sections(f)
end
fig_gcm()
fig_normal_mode()
cs_conv = conv_series()
fig_lmax(cs_conv)

const T_TAB = [0.1, 10.0, 100.0]           # epochs quoted in the table [kyr]

open(joinpath(tables_dir(), "blockB.md"), "w") do io
    println(io, "# Block B — the lateral-viscosity benchmark\n")
    println(io, "Generated by `analysis/run_blockB.jl` from `data/published/bench_visc3d/`. ",
                "Tests A-i, B-i, C-i of the Klemann et al. setup note, ℓmax ", lmax_of(cases),
                ", toroidal coupling on, `visc3d_tol` 1e-3, to 100 kyr as five 20-kyr ",
                "restart segments. There is no reference solution; these are the numbers ",
                "a comparison between codes will be read on.\n")
    for f in BENCH_FORCINGS
        println(io, "## ", uppercasefirst(FORCING_LABEL[f]), "\n")
        println(io, "u_r at the load centre, the peripheral bulge (maximum of u_r along the ",
                    "section), and the geocentre offset |u_CF − u_CM|.\n")
        println(io, "| test | ", join(["u_r at load, $(t) kyr [m]" for t in T_TAB], " | "),
                    " | max u_r at 100 kyr [m] | \\|u_gcm\\| at 100 kyr [m] |")
        println(io, "|---|", repeat("---:|", length(T_TAB) + 2))
        for t in BENCH_TESTS
            haskey(cases, (t, f)) || continue
            c = cases[(t, f)]
            ur = [at_load(c, 5, tk) for tk in T_TAB]
            umax = maximum(c.disp[epoch_rows(c, 100.0), 5])
            gcm = sqrt(sum(abs2, c.gcm[end, 2:4]))
            @printf(io, "| %s-i | %s | %.2f | %.3f |\n", t,
                    join([@sprintf("%.2f", u) for u in ur], " | "), umax, gcm)
            record!("blockB", "$(t)_$(f)_ur_load_100kyr_m" => round(ur[end], digits = 2),
                              "$(t)_$(f)_bulge_100kyr_m"   => round(umax, digits = 2),
                              "$(t)_$(f)_gcm_100kyr_m"     => round(gcm, digits = 3))
        end
        if haskey(cases, ("A", f)) && haskey(cases, ("B", f))
            ra = at_load(cases[("A", f)], 5, 100.0)
            rb = at_load(cases[("B", f)], 5, 100.0)
            @printf(io, "\nB-i / A-i, u_r at the load centre at 100 kyr: %.3f\n", rb / ra)
            record!("blockB", "BA_ratio_$(f)_100kyr" => round(rb / ra, digits = 3))
        end
        println(io)
    end
end
open(joinpath(tables_dir(), "blockB.md"), "a") do io
    println(io, "## A-i against the normal-mode solution\n")
    println(io, "VILMA2 minus the normal-mode solution of M3-L70-V01 at the same truncation ",
                "(ℓmax ", lmax_of(cases), "); maxima along the 0–180° section.\n")
    println(io, "| load | t [kyr] | u_r at load, VILMA2 [m] | normal modes [m] | max\\|Δu_r\\| [m] | max\\|Δu_ϑ\\| [m] | max\\|Δδφ\\| [m² s⁻²] | Δu_gcm [mm] |")
    println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
    for f in BENCH_FORCINGS
        haskey(cases, ("A", f)) || continue
        c = cases[("A", f)]
        for tk in (0.1, 1.0, 10.0, 100.0)
            e = nm_error(c, tk)
            dg = 1e3 * (c.gcm[c.gcm[:, 1] .== tk, 4][1] - gcm_ref(love, σcap, tk, f))
            @printf(io, "| %s | %g | %.3f | %.3f | %.3f | %.3f | %.3f | %+.2f |\n",
                    f, tk, e.ur0, e.ur_nm0, e.ur, e.ut, e.dp, dg)
        end
        es = [nm_error(c, tk) for tk in BENCH_EPOCHS_KYR]
        record!("blockB", "A_$(f)_nm_max_dur_m" => round(maximum(x.ur for x in es), digits = 3),
                          "A_nm_signal_m"       => round(maximum(x.ur_scale for x in es), digits = 1))
    end
    if !isempty(cs_conv)
        println(io, "\n## Truncation\n")
        println(io, "rms Δu_r along the section [m]: A-i against the normal-mode solution to ",
                    "degree ", love.nmax, ", B-i against its own ℓmax 256 run.\n")
        println(io, "| test | load | ℓmax | 1 kyr | 10 kyr | 100 kyr |")
        println(io, "|---|---|---:|---:|---:|---:|")
        ix = [findfirst(==(tk), BENCH_EPOCHS_KYR) for tk in (1.0, 10.0, 100.0)]
        for t in ("A", "B"), f in BENCH_FORCINGS, L in (64, 128, 256)
            haskey(cs_conv, (t, f, L)) || continue
            v = cs_conv[(t, f, L)][ix]
            @printf(io, "| %s-i | %s | %d | %.4f | %.4f | %.4f |\n", t, f, L, v...)
            record!("blockB", "$(t)_$(f)_l$(L)_rms_dur_100kyr_m" => round(v[end], digits = 4))
        end
    end
end
open(joinpath(tables_dir(), "blockB.md"), "a") do io
    rows = [(t, f, bench_cost(t, f)) for t in BENCH_TESTS for f in BENCH_FORCINGS]
    filter!(r -> r[3] !== nothing, rows)
    if !isempty(rows)
        println(io, "\n## Cost\n")
        println(io, "Wall time to 100 kyr at ℓmax ", lmax_of(cases), ", toroidal coupling on, ",
                    "summed over the chain of segment jobs (shared albedo nodes).\n")
        println(io, "| test | load | segments | threads | wall to 100 kyr [h] |")
        println(io, "|---|---|---:|---:|---:|")
        for (t, f, c) in rows
            @printf(io, "| %s-i | %s | %d | %d | %.2f |\n", t, f, c.nseg, c.threads, c.wall_s / 3600)
            record!("blockB", "$(t)_$(f)_wall_100kyr_h" => round(c.wall_s / 3600, digits = 2))
        end
    end
end
record!("blockB", "lmax" => lmax_of(cases), "ncases" => length(cases))
@info "wrote validation/tables/blockB.md"
