# earth_lid.sh — the lid rule of the deglaciation Earth, in one place. Sourced by
# run_deglac.sh (VILMA2's visc3d_lid_depth / visc3d_lid_log10max) and run.sh
# (the same rule written into VILMA1's on-grid viscosity files).
#
# Maxwell elements wholly above LID_DEPTH [m] whose log10 eta exceeds
# LID_LOG10MAX are set to LID_CAP (the clamp ceiling, effectively elastic).
# LID_DEPTH = 0 is off, which is the setting of the reference results: the slow
# part of the held-load relaxation is the cratonic keels at 120-270 km, below
# any lid, so the LGM reference is a fixed 30 kyr hold instead.
LID_DEPTH=0.0
LID_LOG10MAX=22.0
LID_CAP=30.0
