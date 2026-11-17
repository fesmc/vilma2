# Machine configuration: pik_hpc2024 (intended compiler: ifx).
#
# PROJECT-TIER fragment: VILMA2's own copy, which takes precedence over
# configme's shipped machines/pik_hpc2024.mk (orchestrator > user > shipped).
# It reproduces the shipped fragment, with value-safe optimization flags in
# place of its -Ofast, as on albedo (.configme/machines/awi_albedo.mk).
#
# netCDF is auto-detected by configme (load the netCDF module first so
# nf-config/nc-config are on PATH). To pin it, assign INC_NC / LIB_NC here.

# --- Optimization (overrides the shipped -Ofast) -----------------------------
DFLAGS_NODEBUG = -O2 -fp-model precise -march=core-avx2 -mtune=core-avx2 -traceback
