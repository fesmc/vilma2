"""
    Budget

The deglacial sea-level budget: does the water the model delivers equal the
water the ice released?

WHY THIS IS NOT CIRCULAR. The SLE conserves whatever melt source it is handed,
to machine precision, so the run's own `mass_resid` cannot say whether that
source is *right*. Everything here is rebuilt from the run's OUTPUT fields
(`h_ice_eq`, `z_bed_eq`, `h_ice`, `C_ocean`) and compared against a third
quantity the SLE never sees: the above-flotation volume of the ice actually
lost.

    as coded    -(ρi/ρw) ⟨(I − I₀)(1 − C)⟩
    as required -(ρi/ρw) ⟨I(1 − C) − I₀(1 − C₀)⟩
    difference  -(ρi/ρw) ⟨I₀(C − C₀)⟩

The two are identical when the reference is ice-free (I₀ = 0) — which is every
benchmark in the model's suite — and diverge exactly on a cell that carries
grounded marine reference ice, loses it, and floods. That difference was the
marine-reference mask defect, and this is the
measurement that found it. It stays in the pipeline as a standing check, not as
a historical note: it is the only diagnostic here that can catch a melt source
that is wrong rather than merely inconsistent.

Constants are the model's own (`vilma_constants.f90`), not oceanographic
best-estimates. Using 910/1028 instead moves the above-flotation target by about
4 m, which is how sessions 26–27 came to quote a target they later had to
correct.

SIDE OF THE PIPELINE. This module reads `runs/` and is therefore part of the
REDUCTION half: it is included by `analysis/make_reduced.jl` and by nothing
else. The figures read `data/published/` through `analysis/Published.jl`, which
serves the same interface from the reduced files. If you are writing a figure
and reach for this module, you want Published instead.
"""
module Budget

using NCDatasets

export budget_series, budget_at

const RHO_ICE   = 931.0        # vilma_constants.f90:27
const RHO_WATER = 1000.0       # vilma_constants.f90:28
const RR        = RHO_ICE / RHO_WATER

"""Quadrature weights in sin(lat), normalised to 1; handles either ordering."""
function lat_weights(lat)
    n = length(lat); s = sind.(Float64.(lat)); e = zeros(n + 1)
    descending = s[1] > s[end]
    e[1]     = descending ?  1.0 : -1.0
    e[n + 1] = descending ? -1.0 :  1.0
    for j in 2:n
        e[j] = 0.5 * (s[j - 1] + s[j])
    end
    w = [abs(e[j + 1] - e[j]) for j in 1:n] ./ 2
    @assert abs(sum(w) - 1) < 1e-12 "lat weights sum to $(sum(w))"
    w
end

"""
`vilma_sle.f90`'s `ocean_function`: ocean where the bed is below the sea surface
AND the ice there floats rather than grounds. Rebuilt here rather than read from
the file so that C₀ — the reference coastline, which the output does not carry —
is constructed by exactly the expression the model uses.
"""
ocean_function(topo, ice) =
    [(topo[i, j] < 0 && RHO_ICE * ice[i, j] < -RHO_WATER * topo[i, j]) ? 1.0 : 0.0
     for i in axes(topo, 1), j in axes(topo, 2)]

"""
    budget_series(path; nrow = 24) -> NamedTuple

Walk `nrow` time levels of `<run>/out.nc` and return, all in barystatic metres
(sphere mean divided by that moment's ocean fraction):

    time      yr, negative = past
    rsl_mean  ⟨C·rsl⟩/⟨C⟩, the sea level the run actually produced
    coded     the melt source as `vilma_sle.f90` computes it
    required  the same source with each endpoint masked by its own coastline
    target    the above-flotation volume of the reference ice that was lost
    flooding  ⟨z_bed_eq·(C − C₀)⟩, the subgrid term that charges the ocean for
              water filling a newly flooded cell
    ocean_frac

On a model built after the mask fix it is `required + flooding` that reproduces
`rsl_mean`, to ~1e-3 m — that equality is the check that this decomposition
describes the code which actually ran rather than an idealisation of it. Before
the fix it was `coded + flooding` that matched, and the gap between the two
columns was the missing water.
"""
function budget_series(path::AbstractString; nrow::Int = 24)
    ds  = NCDataset(path)
    try
        lat = ds["lat"][:]; lon = ds["lon"][:]; t = ds["time"][:]; nt = length(t)
        heq = ds["h_ice_eq"][:, :]; beq = ds["z_bed_eq"][:, :]
        w = lat_weights(lat); nlon = length(lon)
        sph(F) = sum((sum(@view F[:, j]) / nlon) * w[j] for j in eachindex(w))

        C0  = ocean_function(beq, heq)
        haf = max.(heq .- (RHO_WATER / RHO_ICE) .* max.(-beq, 0.0), 0.0)
        ref_frac = sph(C0)

        rows = unique(round.(Int, range(1, nt, length = min(nrow, nt))))
        n = length(rows)
        out = (time = zeros(n), rsl_mean = zeros(n), coded = zeros(n),
               required = zeros(n), target = zeros(n), flooding = zeros(n),
               ocean_frac = zeros(n), ref_ocean_frac = ref_frac)
        for (k, it) in enumerate(rows)
            h   = ds["h_ice"][:, :, it]
            C   = ds["C_ocean"][:, :, it]
            rsl = ds["rsl"][:, :, it]
            oc  = sph(C)
            lost_af = min.(max.(heq .- h, 0.0), haf)
            out.time[k]       = t[it]
            out.ocean_frac[k] = oc
            out.rsl_mean[k]   = sph(C .* rsl) / oc
            out.coded[k]      = -RR * sph((h .- heq) .* (1 .- C)) / oc
            out.required[k]   = -RR * (sph(h .* (1 .- C)) - sph(heq .* (1 .- C0))) / oc
            out.target[k]     = RR * sph(lost_af) / oc
            out.flooding[k]   = sph(beq .* (C .- C0)) / oc
        end
        out
    finally
        close(ds)
    end
end

"""
    budget_at(path) -> NamedTuple

The present-day row alone — the one the paper quotes. `delivered` is the run's
own ocean-mean rsl; `target` the above-flotation volume lost; `closure` their
ratio. `mask_term` is `coded − required`, the water the endpoint-masking defect
used to discard: the fixed model delivers `required + flooding`, so this is the
size of what the fix recovered, not a residual error.
"""
function budget_at(path::AbstractString)
    b = budget_series(path; nrow = 2)
    k = length(b.time)
    (; time = b.time[k], delivered = b.rsl_mean[k], target = b.target[k],
       closure = b.rsl_mean[k] / b.target[k],
       mask_term = b.coded[k] - b.required[k],
       required_plus_flooding = b.required[k] + b.flooding[k],
       coded_plus_flooding = b.coded[k] + b.flooding[k],
       ocean_frac = b.ocean_frac[k])
end

end # module
