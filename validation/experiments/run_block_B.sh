#!/usr/bin/env bash
#
# run_block_B.sh — block B, the lateral-viscosity benchmark (Klemann et al.
# setup note; docs/benchmarks/visc3d.qmd, driver tests/bench_visc3d.f90).
#
# Every case runs to 100 kyr as a chain of 20-kyr segments, one SLURM job each,
# joined by the driver's restart: segment k starts only when k-1 succeeded
# (afterok), writes restart_<t>kyr.nc, and appends to the protocol files. The
# chain reproduces one uninterrupted run byte for byte. A-i and B-i are
# axisymmetric and take minutes; C-i at lmax 128 with the toroidal field on took
# about 3.7 h per segment on a levante node at 64 threads.
#
#   ./experiments/run_block_B.sh                   # A, B, C x heav, ramp at lmax 128
#   CASES="C" FORCINGS="heav" ./experiments/run_block_B.sh
#   SUBMIT=0 ./experiments/run_block_B.sh          # stage only
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib_stage.sh
. "$HERE/lib_stage.sh"
RUNS="${RUNS:-$VALIDATION/runs/blockB}"
STAGE="$VALIDATION/experiments/.stage"
TMPL="$HERE/nml/bench_visc3d.nml.tmpl"

SUBMIT=${SUBMIT:-1}
LMAX=${LMAX:-128}
LTOR=${LTOR:-.true.}
V3DTOL=${V3DTOL:-1.0e-3}
CASES=${CASES:-"A B C"}
FORCINGS=${FORCINGS:-"heav ramp"}
SEG_YR=${SEG_YR:-20000}
T_END_YR=${T_END_YR:-100000}
OMP=${OMP:-$V_OMP_DEGLAC}

mkdir -p "$RUNS" "$STAGE"

# The driver runs from $MODEL/bin via the runme alias `visc3d`.
( cd "$MODEL" && make bench_visc3d openmp=1 >/dev/null )
[ -x "$MODEL/bin/bench_visc3d.x" ] || { echo "run_block_B: bench_visc3d.x was not built" >&2; exit 1; }

echo "model   : $MODEL ($(cd "$MODEL" && git rev-parse --short HEAD))    machine: $V_MACHINE"
echo "runs  : $RUNS    lmax=$LMAX  l_toroidal=$LTOR  segments of $SEG_YR yr to $T_END_YR yr"
echo

# nml <out> <test> <forcing> <t_end> <restart_in> <restart_out>
nml() {
  sed -e "s/@TEST@/$2/g" -e "s/@FORCING@/$3/g" -e "s/@TEND@/$4/g" \
      -e "s|@RIN@|$5|g" -e "s|@ROUT@|$6|g" -e "s/@LMAX@/$LMAX/g" \
      -e "s/@LTOR@/$LTOR/g" -e "s/@V3DTOL@/$V3DTOL/g" "$TMPL" > "$1"
  if grep -qE '@[A-Z0-9_]+@' "$1"; then
    echo "run_block_B: unsubstituted placeholder in $1" >&2; exit 1
  fi
}

tag_tor=""; [ "$LTOR" = ".false." ] && tag_tor=".notor"
nseg=$(( T_END_YR / SEG_YR ))

for T in $CASES; do
  # C is fully 3-D; A and B are axisymmetric (mmax = 0) and take minutes.
  wall=01:00:00; [ "$T" = "C" ] && wall=06:00:00
  for F in $FORCINGS; do
    id="${T}-i_${F}.l${LMAX}${tag_tor}"
    dir="$RUNS/$id"
    echo ">>> $id   ($nseg segments, wall=$wall each)"

    # Stage once with segment 1, so runme links the binary, input and data.
    nml "$STAGE/bv3d_${id}_seg1.nml" "$T" "$F" "${SEG_YR}.0" "" "restart_$((SEG_YR/1000))kyr.nc"
    v_stage "$dir" visc3d "$STAGE/bv3d_${id}_seg1.nml" "$wall" "$OMP" spread

    # One namelist and one job script per segment, in the run directory.
    seg1=$(basename "$STAGE/bv3d_${id}_seg1.nml")
    prev=""
    for k in $(seq 1 "$nseg"); do
      t1=$(( k*SEG_YR )); rin=""; rout="restart_$((t1/1000))kyr.nc"
      [ "$k" -gt 1 ] && rin="restart_$(( (k-1)*SEG_YR/1000 ))kyr.nc"
      [ "$k" -eq "$nseg" ] && rout=""
      nml "$dir/seg$k.nml" "$T" "$F" "${t1}.0" "$rin" "$rout"
      sed -e "s|^srun \(.*\) $seg1|srun \1 seg$k.nml|" \
          -e "s|^#SBATCH --output=.*|#SBATCH --output=./out.seg$k.out|" \
          -e "s|^#SBATCH --error=.*|#SBATCH --error=./out.seg$k.err|" \
          -e "s|^#SBATCH --job-name=.*|#SBATCH --job-name=B${T}${F:0:1}$k|" \
          "$dir/job.submit" > "$dir/job.seg$k.submit"
      grep -q "seg$k.nml" "$dir/job.seg$k.submit" || { echo "run_block_B: srun line not rewritten in $dir/job.seg$k.submit" >&2; exit 1; }
      if [ "$SUBMIT" = "1" ]; then
        dep=""; [ -n "$prev" ] && dep="--dependency=afterok:$prev"
        prev=$(cd "$dir" && sbatch --parsable $dep "job.seg$k.submit")
        echo "    seg$k -> t=${t1} yr  job $prev"
      fi
    done
  done
done

echo
echo "done. squeue -u $USER"
