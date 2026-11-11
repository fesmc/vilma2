"""
    NormalMode

The exact solution of benchmark test A-i (the √-profile cap on a laterally
uniform M3-L70-V01) from normal-mode loading Love numbers — the reference the
lateral-viscosity benchmark itself does not provide.

For an incompressible layered Maxwell sphere the Heaviside response of each
degree is a finite sum of decaying exponentials. The table
(`data/published/bench_visc3d/mod_M3-L70-V01`, vendored in the model repo from
giapy, Spada et al. 2011; model `data/benchmarks/PROVENANCE.md`) gives, per
degree n = 1…256, the elastic numbers (k_e, h_e, l_e) and nine modes (s_i [1/kyr],
amplitudes a_k, a_h, a_l) with

    x(t) = x_e + Σ_i a_i (1 − e^{s_i t}),     x_f = x_e + Σ_i a_i,

which the reader checks against the table's own fluid line. Degree 1 is in the
centre-of-mass frame (k_1 = −1), the frame the benchmark driver runs in. Degree
0 has h = l = k = 0 (incompressible): the load's own potential only.

Fields at the surface, for a zonal load σ = Σ σ_n P_n(cos ϑ):

    Φ_n/g = 3 σ_n / (ρ̄ (2n+1))
    u_r = Σ h_n Φ_n/g P_n,   u_ϑ = Σ l_n Φ_n/g ∂_ϑ P_n,   δφ = g N = Σ (f+k_n) Φ_n P_n

where f(t) is the fraction of the load in place (1, or t/10 kyr on the ramp):
the load's own potential follows the load, not a Love-number convolution.
    u_CF − u_CM (along the load axis) = (h_1 + 2 l_1)/3 · Φ_1/g

with ρ̄ and g from the table's own earth model and the model's G.
"""
module NormalMode

using Printf

export LoveTable, read_love_table, cap_legendre, love_eff, section_ref, gcm_ref

const G_NEWTON = 6.67430e-11        # the model's grav_G (vilma_constants.f90)

struct LoveTable
    nmax::Int
    a::Float64                      # surface radius [m]
    rhobar::Float64                 # mean density [kg m^-3]
    g0::Float64                     # surface gravity [m s^-2]
    e::Matrix{Float64}              # (3, nmax): k_e, h_e, l_e for n = 1..nmax
    s::Matrix{Float64}              # (nmodes, nmax) [1/kyr]
    amp::Array{Float64,3}           # (3, nmodes, nmax): a_k, a_h, a_l
end

"Parse the giapy/TABOO table and check each degree's modes sum to its fluid line."
function read_love_table(path)
    lines = [split(strip(l)) for l in eachline(path) if !isempty(strip(l))]
    lay = [parse.(Float64, lines[i]) for i in 1:5]          # layer r rho mu eta
    r   = [l[2] for l in lay];  rho = [l[3] for l in lay]
    rb  = vcat(r[2:end], 0.0)
    M   = sum(4π / 3 * rho[i] * (r[i]^3 - rb[i]^3) for i in eachindex(r))
    a   = r[1]
    E = Vector{Vector{Float64}}(); S = Vector{Vector{Float64}}(); A = Vector{Matrix{Float64}}()
    i = 6
    while i <= length(lines)
        h = parse.(Float64, lines[i]);  n, nm = Int(h[1]), Int(h[2])
        n == length(E) + 1 || error("$path: degree $n out of order")
        md = [parse.(Float64, lines[i + j]) for j in 1:nm]
        f  = parse.(Float64, lines[i + nm + 1])
        amp = reduce(hcat, [m[3:5] for m in md])              # (3, nm)
        fl  = h[3:5] .+ vec(sum(amp, dims = 2))
        maximum(abs.(fl .- f[3:5]) ./ max.(abs.(f[3:5]), 1e-12)) < 1e-3 ||
            error("$path: degree $n modes do not sum to the fluid line ($fl vs $(f[3:5]))")
        push!(E, h[3:5]);  push!(S, [m[2] for m in md]);  push!(A, amp)
        i += nm + 2
    end
    nmax = length(E);  nm = length(S[1])
    amp = zeros(3, nm, nmax)
    for n in 1:nmax
        amp[:, :, n] = A[n]
    end
    LoveTable(nmax, a, 3M / (4π * a^3), G_NEWTON * M / a^2,
              reduce(hcat, E), reduce(hcat, S), amp)
end

"n-point Gauss–Legendre nodes and weights on [0, 1] (Newton on P_n)."
function gauss_legendre01(n)
    x = zeros(n); w = zeros(n)
    for i in 1:(n + 1) ÷ 2
        z = cos(π * (i - 0.25) / (n + 0.5)); pp = 0.0
        for _ in 1:100
            p1, p2 = 1.0, 0.0
            for j in 1:n
                p1, p2 = ((2j - 1) * z * p1 - (j - 1) * p2) / j, p1
            end
            pp = n * (z * p1 - p2) / (z^2 - 1)
            z1 = z; z -= p1 / pp
            abs(z - z1) < 1e-15 && break
        end
        x[i] = 0.5(1 - z); x[n + 1 - i] = 0.5(1 + z)
        w[i] = 1 / ((1 - z^2) * pp^2); w[n + 1 - i] = w[i]
    end
    x, w
end

"""
    cap_legendre(nmax; alpha_deg = 10, h = 1500, rho = 931) -> σ_n, n = 0..nmax

Legendre coefficients [kg m⁻²] of the √-profile cap σ(ϑ) = ρ h √((cos ϑ − c)/(1 − c)),
c = cos α, by the driver's own quadrature: x = c + (1 − c)s² makes the √ edge a
polynomial, exact for nq ≥ nmax + 2.
"""
function cap_legendre(nmax; alpha_deg = 10.0, h = 1500.0, rho = 931.0)
    c = cosd(alpha_deg); s, w = gauss_legendre01(nmax + 4)
    σ = zeros(nmax + 1); p = zeros(nmax + 1)
    for q in eachindex(s)
        x = c + (1 - c) * s[q]^2
        legendre!(p, x)
        σ .+= w[q] * s[q] * 2(1 - c) * s[q] .* p
    end
    [(2n + 1) / 2 * rho * h * σ[n + 1] for n in 0:nmax]
end

"""
    love_eff(tab, n, t_kyr, forcing; t_ramp = 10) -> (k, h, l)

Effective Love numbers at time t: the Heaviside response for "heav", its
convolution with a linear ramp over 0–t_ramp kyr for "ramp". Degree 0 is zero.
"""
function love_eff(tab::LoveTable, n, t, forcing; t_ramp = 10.0)
    n == 0 && return (0.0, 0.0, 0.0)
    xe = tab.e[:, n]; s = tab.s[:, n]; A = tab.amp[:, :, n]
    if forcing == "heav"
        x = xe .+ A * (1 .- exp.(s .* t))
    elseif forcing == "ramp"
        m  = min(t, t_ramp)
        xf = xe .+ vec(sum(A, dims = 2))
        x  = (m .* xf .- A * ((exp.(s .* t) .- exp.(s .* (t - m))) ./ s)) ./ t_ramp
    else
        error("forcing must be heav or ramp")
    end
    (x[1], x[2], x[3])
end

"""
    section_ref(tab, σ, dist_deg, t_kyr, forcing; nmax) -> (u_r, u_ϑ, δφ)

The fields at angular distances `dist_deg` from the load centre, summed to degree
`nmax` (the run's ℓmax, or the table's 256).
"""
function section_ref(tab::LoveTable, σ, dist, t, forcing; nmax = tab.nmax)
    nmax <= min(tab.nmax, length(σ) - 1) || error("nmax $nmax beyond the table or the load")
    L = [love_eff(tab, n, t, forcing) for n in 0:nmax]
    φg = [3σ[n + 1] / (tab.rhobar * (2n + 1)) for n in 0:nmax]    # Φ_n / g
    ch = [φg[n + 1] * L[n + 1][2] for n in 0:nmax]
    cl = [φg[n + 1] * L[n + 1][3] for n in 0:nmax]
    f  = load_fraction(t, forcing)
    ck = [φg[n + 1] * (f + L[n + 1][1]) * tab.g0 for n in 0:nmax]
    ur = zeros(length(dist)); ut = zeros(length(dist)); dp = zeros(length(dist))
    p = zeros(nmax + 1); q = zeros(nmax + 1)
    for (j, d) in enumerate(dist)
        legendre!(p, cosd(d)); legendre1!(q, cosd(d), sind(d))
        ur[j] = sum(ch .* p); ut[j] = sum(cl .* q); dp[j] = sum(ck .* p)
    end
    ur, ut, dp
end

"Fraction of the load in place at t [kyr]."
load_fraction(t, forcing; t_ramp = 10.0) = forcing == "heav" ? 1.0 : min(t / t_ramp, 1.0)

"P_n(x), n = 0..length(p)-1, in place."
function legendre!(p, x)
    p[1] = 1.0
    length(p) > 1 && (p[2] = x)
    for n in 2:length(p) - 1
        p[n + 1] = ((2n - 1) * x * p[n] - (n - 1) * p[n - 1]) / n
    end
    p
end

"∂_ϑ P_n = P_n^1 (Condon–Shortley phase), n = 0..length(q)-1, at x = cos ϑ, sn = sin ϑ."
function legendre1!(q, x, sn)
    q[1] = 0.0
    length(q) > 1 && (q[2] = -sn)
    for n in 2:length(q) - 1
        q[n + 1] = ((2n - 1) * x * q[n] - n * q[n - 1]) / (n - 1)
    end
    q
end

"Geocentre offset u_CF − u_CM along the load axis [m]."
function gcm_ref(tab::LoveTable, σ, t, forcing)
    k, h, l = love_eff(tab, 1, t, forcing)
    (h + 2l) / 3 * σ[2] / tab.rhobar          # Φ_1/g = 3σ_1/(3ρ̄)
end

end # module
