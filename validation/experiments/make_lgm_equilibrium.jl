#!/usr/bin/env julia
#
# make_lgm_equilibrium.jl — VILMA2's OWN LGM equilibrium as the i_eq=2
# reference state, read from a RESTART FILE.
#
#   julia --project=../analysis experiments/make_lgm_equilibrium.jl [run_dir | restart.nc]
#
# WHY THIS REPLACES make_lgm_reference.jl.
#
# `make_lgm_reference.jl` handed the model GLAC-1D's own LGM bed and LGM ice and
# declared them an equilibrium. They are not one -- that bed came from a
# different Earth model -- and i_eq=2 has no way to find out, because the
# reference state is an ASSERTION, not a computation. Worse, the spin-up cannot
# repair it: solid_earth_spinup relaxes the departure of the load from the i_eq
# reference, and with i_eq=2 the reference IS the start-slice state, so the
# departure is zero and the spin-up converges on pass 1 with
# `d|z_bed| = 0.0000 m` and does nothing. Measured:
# VILMA2 then returns 43 % of the deglacial barystatic budget and its
# ocean-mean sea level peaks at 10 ka and falls 12 m while the ice is constant.
#
# The equilibrium VILMA2 needs is one IT computes, under the LGM load, with
# its own Earth structure. Initialization I1 does exactly that: hold the 26 ka
# load against the RTopo-2 datum (equil_time_max = 100 kyr, equil_rate_tol =
# 1e-3, pre_spinup_1d = .true.) until the bed stops moving. That IS the whole
# experiment -- a constant boundary condition relaxed to a tolerance.
#
# WHY A RESTART FILE, NOT out.nc. This script used to read
# the t = -26000 slice of the source run's `out.nc`, which meant the reference
# run had to be a full -26000 -> 0 deglaciation whose 260 remaining steps were
# computed and thrown away -- about 85 % of its wall time, on the one run the
# whole pipeline serialises behind. The state it took is exactly the state the
# spin-up already wrote to `spinup/vilma_restart.nc`: verified bit-identical on
# z_bed, h_ice, rsl and C_ocean, and on bsl to all digits. A restart file is also
# the expected place to take a model state from, so the reference can now be
# rebuilt from any restart -- a spin-up, the end of a run, or a rung of a cascade
# -- rather than from one particular transient's output.
#
# NOTE the `time` in a spin-up restart is the spin-up's OWN relaxation clock (the
# run this was written from converged after 14503 yr of held load), not a forcing
# epoch. It says how long the relaxation took, and nothing about the LGM.
#
# THE DATUM, which is the trap that started all of this. vilma_sle reads topo0 as
# bed elevation relative to the CONTEMPORANEOUS sea surface (vilma_sle.f90:200,
# "bathymetry zeta(0) = topo0, rsl = 0 at t0"). In the source run
# z_bed = z_bed_eq - rsl with rsl a RELATIVE sea level, so z_bed is already the
# bed relative to the LGM sea surface -- no offset is applied here, and none
# should be. (The measured lowering in that state is -100.2 m of ocean-mean rsl,
# against the 104.6 m the ice budget implies.)
#
# The output is written on the model's own Gauss grid, which read_ref2d takes
# directly when its lon/lat dims match nphi/nlat -- so no remap, and the
# reference is bit-exactly the state the source run was in.
#
# Viscosity matters for HOW LONG the spin-up takes, not for the relaxed state
# itself (a Maxwell body forgets its viscosity at full relaxation), but the cap
# is 100 kyr, so use a source run whose Earth structure matches the transient
# you intend to run from this reference.

using NCDatasets
using Printf

const REPO = dirname(@__DIR__)
const ARG  = length(ARGS) >= 1 ? ARGS[1] :
             joinpath(REPO, "runs", "blockD_ref", "D-V2-3D.l170.cold.ref")
# Optional second argument: where to write, for a diagnostic export that must
# not replace the production reference.
const OUT  = length(ARGS) >= 2 ? ARGS[2] :
             joinpath(REPO, "data", "reference", "vilma2_lgm_eq.nc")

# Accept a restart file directly, or a run directory -- in which case prefer the
# spin-up restart (the relaxed state this script wants) over the end-of-run one.
function resolve_restart(a)
    isfile(a) && return a
    isdir(a) || error("$a is neither a restart file nor a directory")
    for c in (joinpath(a, "spinup", "vilma_restart.nc"), joinpath(a, "final", "vilma_restart.nc"))
        isfile(c) && return c
    end
    error("no restart under $a: expected spinup/vilma_restart.nc or final/vilma_restart.nc.\n" *
          "A reference run writes the spin-up restart when the held load reaches\n" *
          "equil_rate_tol; if it is missing the spin-up did not converge.")
end

const SRC = resolve_restart(ARG)
ds = NCDataset(SRC)
haskey(ds, "z_bed") || error("$SRC has no z_bed: this is not a coupled-state restart")
nt = length(ds["time"][:])
nt == 1 || @warn "restart carries $nt time records; taking the last"
k  = nt
@info @sprintf("source %s (relaxation clock %.1f yr)", SRC, Float64(ds["time"][k]))

lon = Float64.(ds["lon"][:]); lat = Float64.(ds["lat"][:])
bed = Float64.(ds["z_bed"][:, :, k])      # equilibrium bed, vs the LGM sea surface
ice = Float64.(ds["h_ice"][:, :, k])      # the LGM load, as the model remapped it
rsl = Float64.(ds["rsl"][:, :, k])
C   = Float64.(ds["C_ocean"][:, :, k])
bsl = Float64(ds["bsl"][k])
close(ds)

w  = repeat(reshape(cosd.(lat), 1, :), length(lon), 1)
om = sum(rsl .* C .* w) / sum(C .* w)
@printf("equilibrium state: ocean-mean rsl %.2f m, bsl %.2f m, ocean fraction %.4f\n",
        om, bsl, sum(C .* w) / sum(w))
@printf("bed  [%.1f, %.1f] m\nice  [%.1f, %.1f] m\n", extrema(bed)..., extrema(ice)...)

mkpath(dirname(OUT)); isfile(OUT) && rm(OUT)
d = NCDataset(OUT, "c")
defDim(d, "lon", length(lon)); defDim(d, "lat", length(lat))
defVar(d, "lon", lon, ("lon",), attrib = ["units" => "degrees_east",  "axis" => "X"])
defVar(d, "lat", lat, ("lat",), attrib = ["units" => "degrees_north", "axis" => "Y"])
defVar(d, "bedrock_topography", bed, ("lon", "lat"), attrib = [
    "long_name" => "VILMA2's own equilibrium bed under the LGM load, " *
                   "relative to the sea surface at that epoch", "units" => "m"])
defVar(d, "ice_thickness", ice, ("lon", "lat"), attrib = [
    "long_name" => "grounded-ice thickness under the LGM load, as remapped by the model",
    "units" => "m"])
d.attrib["source_restart"] = SRC
d.attrib["comment"] = "i_eq=2 reference computed BY VILMA2 (initialization I1 " *
    "against the RTopo-2 datum), not asserted from an external bed, and read from " *
    "the spin-up restart rather than a transient slice. See " *
    "validation/experiments/make_lgm_equilibrium.jl."
close(d)
@info "wrote $OUT"
