# Machine configuration: dkrz_levante (intended compiler: ifx).
#
# PROJECT-TIER fragment: VILMA2's own copy, which takes precedence over
# configme's shipped machines/dkrz_levante.mk (orchestrator > user > shipped).
# It reproduces the shipped fragment, with value-safe optimization flags in
# place of its -Ofast, as on albedo (.configme/machines/awi_albedo.mk).
#
# netCDF is auto-detected by configme (load the netCDF module first so
# nf-config/nc-config are on PATH). To pin it, assign INC_NC / LIB_NC here.

# --- Optimization (overrides the shipped -Ofast) -----------------------------
# Levante compute nodes are AMD EPYC 7763 (Milan/Zen3): AVX2, no AVX-512.
# (Note: the login-node gcc is too old to name Zen3 and mis-resolves
# -march=native to znver1 -- confirm the CPU with lscpu, not gcc.)
DFLAGS_NODEBUG = -O2 -fp-model precise -march=znver3 -traceback
