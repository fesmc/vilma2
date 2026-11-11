# DKRZ levante. Sourced by experiments/lib_stage.sh.
V_ACCOUNT=ba1442
V_PARTITION=compute
V_NCORE=128                 # cores per node
V_EXCLUSIVE=1               # whole node, all its memory; threads set by env only
V_OMP_DEGLAC=64
V_WALL_DEGLAC=08:00:00
V_VILMA1_ROOT=/work/ba1442/robinson/models/vilma1
V_VILMA1_INPUT=/work/ba1442/robinson/models/climber-x/input/vilma
V_VISC_BAGGE="${V_VISC_BAGGE:-$VALIDATION/data/input/bagge2021.nc}"

v_queue() { echo compute; }
