# AWI albedo. Sourced by experiments/lib_stage.sh.
V_ACCOUNT=envi.p_forclima
V_PARTITION=smp
V_NCORE=128                 # cores per node
V_EXCLUSIVE=0               # 0: allocate only the threads a job uses (shared node)
V_OMP_DEGLAC=64             # threads of the deglaciation and bench_visc3d runs
V_WALL_DEGLAC=08:00:00
V_VILMA1_ROOT="$HOME/models/climber-x-v24/src/vilma1"       # include/, lib/vega_pism.a
V_VILMA1_INPUT="$HOME/models/climber-x-v24/input/vilma"     # static VILMA1 input files
# Bagge et al. (2021) on lon, lat, r, quantised to 0.01 dex: the file the model
# repo vendored as input/bagge2021.nc until 2676879.
V_VISC_BAGGE="${V_VISC_BAGGE:-$VALIDATION/data/input/bagge2021.nc}"

# runme queue alias (= QOS on albedo) for a wall time hh:mm:ss
v_queue() {
  local h=${1%%:*} m=${1#*:}; m=${m%%:*}
  local min=$(( 10#$h * 60 + 10#$m ))
  if   [ "$min" -le 30 ];  then echo 30min
  elif [ "$min" -le 720 ]; then echo 12h
  else echo 48h; fi
}
