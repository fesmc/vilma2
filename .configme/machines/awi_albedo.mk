# Machine configuration: awi_albedo (intended compiler: ifx).
#
# PROJECT-TIER fragment: VILMA2's own copy, which takes precedence over
# configme's shipped machines/awi_albedo.mk (orchestrator > user > shipped). It
# reproduces the shipped fragment, with value-safe optimization flags in place
# of its -Ofast, and refuses ifx 2024.1.
#
# netCDF is auto-detected by configme (load the netCDF module first so
# nf-config/nc-config are on PATH). To pin it, assign INC_NC / LIB_NC here.

# --- Optimization (overrides the shipped -Ofast) -----------------------------
# Albedo compute nodes (mpp/smp) are uniformly AMD EPYC 7702 (Rome/Zen2), AVX2
# with no AVX-512, so znver2 is a pure tuning win with no ISA risk.
# -fp-model precise instead of -Ofast: the coupled kryos stack builds VILMA2 with
# these flags, since -Ofast breaks the round-off conservation of kryos'
# test_globe_ice; VILMA2 builds and runs standalone with the same flags.
DFLAGS_NODEBUG = -O2 -fp-model precise -march=znver2 -traceback

# --- Compiler version --------------------------------------------------------
# Build with intel-oneapi-compilers/2024.0.0. ifx 2024.1 miscompiles (an inlined
# function result passed as an argument loses a WHERE; it broke kryos'
# test_geo_vilma2), and it is the login default on albedo, so the build stops
# on it. clean, usage and showconfig still run.
ifneq ($(filter-out clean usage showconfig,$(or $(MAKECMDGOALS),usage)),)
ifneq ($(findstring (IFX) 2024.1.,$(shell $(FC) --version 2>/dev/null)),)
$(error ifx 2024.1 miscompiles on albedo. Run: module unload intel-oneapi-compilers; module load intel-oneapi-compilers/2024.0.0; make clean)
endif
endif
