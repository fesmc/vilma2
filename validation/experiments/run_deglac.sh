#!/usr/bin/env bash
#
# run_deglac.sh — stage and submit the last-deglaciation runs (GLAC-1D,
# 26 ka to present, PREM geometry with the Bagge et al. 2021 viscosity).
#
#   ./experiments/run_deglac.sh blockD       VILMA2, Earth structures $STRUCTS
#   ./experiments/run_deglac.sh blockD-res   VILMA2 3-D at the truncations $LMAX_LIST_D
#   ./experiments/run_deglac.sh blockD-tor   VILMA2 3-D with the toroidal field (V7)
#   ./experiments/run_deglac.sh blockDV1     the VILMA1 backend, $STRUCTS (needs EXE)
#   ./experiments/run_deglac.sh blockE-tol   VILMA2 3-D over visc_3d_tol
#
# run.sh sets the environment of each stage (reference state, frame, scheme);
# SUBMIT=0 stages without submitting. Templates: nml/blockD.nml.tmpl and
# nml/blockDV1.nml.tmpl.
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib_stage.sh
. "$HERE/lib_stage.sh"
# shellcheck source=earth_lid.sh
. "$HERE/earth_lid.sh"

STAGE="$VALIDATION/experiments/.stage"
TMPL_D="$HERE/nml/blockD.nml.tmpl"
TMPL_DV1="$HERE/nml/blockDV1.nml.tmpl"
SUBMIT=${SUBMIT:-1}
OMP=${OMP:-$V_OMP_DEGLAC}
WALL=${WALL:-$V_WALL_DEGLAC}
mkdir -p "$STAGE"

# check_nml <file>: an unsubstituted placeholder would only fail in the job.
check_nml() {
  if grep -qE '@[A-Z0-9_]+@' "$1"; then
    echo "unsubstituted placeholder(s) in $1:" >&2
    grep -oE '@[A-Z0-9_]+@' "$1" | sort -u | sed 's/^/  /' >&2
    exit 1
  fi
}

# nmlD <name> <lmax> <scheme> <rtol> <equil> <pre1d> <v3d> <v3dtol> <t0> [t1]
# DEG1 defaults to the production value (cm frame).
nmlD() {
  local out="$STAGE/$1.nml"
  sed -e "s/@LMAX@/$2/g"   -e "s/@SCHEME@/$3/g" -e "s/@RTOL@/$4/g" \
      -e "s/@EQUIL@/$5/g"  -e "s/@PRE1D@/$6/g"  -e "s/@V3D@/$7/g" \
      -e "s/@V3DTOL@/$8/g" -e "s/@T0@/$9/g"     -e "s/@T1@/${10:-0.0}/g" \
      -e "s/@EQTOL@/${EQTOL:-1.0e-3}/g" \
      -e "s/@IEQ@/${IEQ:-1}/g" -e "s|@ZEQ@|${ZEQ:-}|g" -e "s|@HEQ@|${HEQ:-}|g" \
      -e "s|@VISCF@|${VISCF:-$V_VISC_BAGGE}|g" \
      -e "s|@DEG1@|${DEG1:-cm}|g" -e "s|@LTOR@|${LTOR:-.false.}|g" \
      -e "s|@LIDDEPTH@|$LID_DEPTH|g" -e "s|@LIDMAX@|$LID_LOG10MAX|g" \
      "$TMPL_D" > "$out"
  check_nml "$out"
  printf '%s' "$out"
}

# nmlDV1 <name> <lmax> <v3d> <t0> [t1]
# The VILMA1 backend has its own template: blockD.nml.tmpl sets no `solver`, so
# a run staged from it with the VILMA1 binary is still a VILMA2 run.
nmlDV1() {
  local out="$STAGE/$1.nml"
  sed -e "s/@LMAX@/$2/g"  -e "s/@V3D@/$3/g" -e "s/@T0@/$4/g" \
      -e "s/@T1@/${5:-0.0}/g" \
      -e "s|@VINDIR@|$VINDIR|g" \
      -e "s|@VISC1D@|${VISC1D:-visko.inp}|g" \
      -e "s|@VISC3D@|$VISC3D|g" \
      -e "s/@NSUB@/${NSUB:-12}/g"   -e "s/@LPREM@/${LPREM:-1}/g" \
      -e "s/@IEQ@/${IEQ:-1}/g" -e "s|@ZEQ@|${ZEQ:-}|g" -e "s|@HEQ@|${HEQ:-}|g" \
      "$TMPL_DV1" > "$out"
  check_nml "$out"
  grep -qE '^\s*solver\s*=\s*"v1"' "$out" || {
    echo "nmlDV1: $out does not set solver=\"v1\"" >&2; exit 1; }
  printf '%s' "$out"
}

# go <id> <nml> [alt-binary]
go() {
  local id=$1 par=$2 alt=${3:-}
  local dir="$RUNS/$id"
  echo ">>> $id   (wall=$WALL, threads=$OMP${alt:+, exe=$(basename "$alt")})"
  v_stage "$dir" bin/vilma.x "$par" "$WALL" "$OMP" spread
  [ -n "$alt" ] && cp "$alt" "$dir/vilma.x"
  if [ "$SUBMIT" = "1" ]; then
    ( cd "$dir" && sbatch --job-name="${id%%.*}" job.submit )
  fi
}

mode=${1:-}
case "$mode" in
  blockD|blockD-res) RUNS="${RUNS_BLOCKD:-$VALIDATION/runs/blockD_v2}" ;;
  blockDV1)          RUNS="${RUNS_BLOCKD:-$VALIDATION/runs/blockD_v1}" ;;
  blockD-tor)        RUNS="${RUNS_BLOCKD:-$VALIDATION/runs/blockD_tor}" ;;
  blockE-tol)        RUNS="${RUNS_BLOCKD:-$VALIDATION/runs/blockE}" ;;
  *) echo "usage: $0 {blockD|blockD-res|blockDV1|blockD-tor|blockE-tol}" >&2; exit 1 ;;
esac
mkdir -p "$RUNS"
echo "model   : $MODEL ($(cd "$MODEL" && git rev-parse --short HEAD))"
echo "machine : $V_MACHINE    mode: $mode    runs: $RUNS"
echo

case "$mode" in

  # VILMA2 at lmax $LMAX_D, forward Euler (cfl 1), by Earth structure:
  #   1Da  l_visc_3d=.false.                  the PREM/"Clemens" radial profile
  #   1Db  l_visc_3d=.true., visc_3d_tol=99    the laterally uniform file
  #        (make_uniform_visc.jl: per-level geometric mean of the clamped Bagge
  #        field), which D-V1-1Db reads too
  #   3D   l_visc_3d=.true., visc_3d_tol=1e-3  full lateral variation
  blockD)
    L=${LMAX_D:-170}; SCH=${SCHEME:-fe}; RT=${RTOL:-1.0e-4}
    UNIF=${UNIF_FE:-$VALIDATION/data/reference/bagge2021_uniform_vilma2.nc}
    for st in ${STRUCTS:-3D 1Db}; do
      case $st in
        1Da) v3=.false.; tol=1.0e-3; pre=.false.; vf=$V_VISC_BAGGE ;;
        1Db) v3=.true.;  tol=99.0;   pre=.true.;  vf=$UNIF ;;
        3D)  v3=.true.;  tol=1.0e-3; pre=.true.;  vf=$V_VISC_BAGGE ;;
        *) echo "unknown structure $st" >&2; exit 1 ;;
      esac
      [ -f "$vf" ] || { echo "blockD: $vf missing" >&2; exit 1; }
      # Exported, not a prefix: nmlD expands in go's argument list.
      export VISCF="$vf"
      go "${MODEL_TAG:-D-V2}-$st.l$L.$SCH" \
         "$(nmlD "${NML_TAG:-dF}_${st}_${L}" "$L" "$SCH" "$RT" "${EQUIL:-100000.0}" \
            "${PRE1D:-$pre}" "$v3" "$tol" -26000.0)" ${EXE:-}
    done
    ;;

  # The VILMA1 backend on the same reference file, ice and output grid as blockD.
  # Matched to the D-V2 rows:
  #   * vilma1_input_vilma2grid/densi.inp with vilma_v1_l_prem = 0: VILMA1 reads
  #     VILMA2's 11-layer table instead of building its polynomial PREM;
  #   * the viscosity from 3-D files resampled onto that mesh's element centres
  #     (make_vilma1_visc_ongrid.jl; VILMA1 maps the levels by index): the
  #     clamped Bagge field for 3D, the laterally uniform field for 1Db.
  blockDV1)
    L=${LMAX_D:-170}
    [ -n "${EXE:-}" ] || { echo "blockDV1: set EXE to a vilma_v1=1 binary" >&2; exit 1; }
    export VINDIR="${VINDIR:-$HERE/vilma1_input_vilma2grid}" LPREM=0 VISC1D=visko_bagge.inp
    for st in ${STRUCTS:-3D 1Db}; do
      case $st in
        1Db) v3=.true.; f=unif_vilma2grid.nc  ;;
        3D)  v3=.true.; f=bagge_vilma2grid.nc ;;
        *) echo "unknown structure $st for blockDV1" >&2; exit 1 ;;
      esac
      [ -f "$VINDIR/$f" ] || { echo "blockDV1: $VINDIR/$f missing -- see make_vilma1_visc_ongrid.jl" >&2; exit 1; }
      python3 "$HERE/check_vilma1_visc_grid.py" "$VINDIR/densi.inp" "$VINDIR/$f"
      export VISC3D="$f"
      go "${MODEL_TAG:-D-V1}-$st.l$L.${SCHEME:-fe}" \
         "$(nmlDV1 "${NML_TAG:-dV1}_${st}_${L}" "$L" "$v3" -26000.0)" "$EXE"
    done
    ;;

  # V7: the D-V2-3D row with the toroidal field of the 3-D elements carried
  # (l_toroidal); everything else identical, so tor minus D-V2-3D.l170.fe is the
  # toroidal effect. Six more memory channels per 3-D element: allow more wall.
  blockD-tor)
    L=${LMAX_D:-170}; SCH=${SCHEME:-fe}; RT=${RTOL:-1.0e-4}
    export VISCF="$V_VISC_BAGGE" LTOR=.true.
    go "${MODEL_TAG:-D-V2}-3D.l$L.tor" \
       "$(nmlD "${NML_TAG:-dT}_tor_${L}" "$L" "$SCH" "$RT" "${EQUIL:-100000.0}" \
          "${PRE1D:-.true.}" .true. 1.0e-3 -26000.0)"
    ;;

  # visc_3d_tol, the speed/accuracy lever of the 3-D path, on the D-V2-3D row.
  # tol = 0 keeps every Maxwell element with any lateral spread 3-D.
  blockE-tol)
    L=${LMAX_E3:-170}
    for tol in ${TOL_LIST:-1.0e-3 1.0e-2 1.0e-1 0.0}; do
      case $tol in
        1.0e-3) tag=1e-3 ;;   # default
        1.0e-2) tag=1e-2 ;;
        1.0e-1) tag=1e-1 ;;
        0.0)    tag=all  ;;
        *)      tag=$(echo "$tol" | tr -d '.') ;;
      esac
      go "E3-tol$tag.l$L" \
         "$(nmlD "dE3_${tag}_${L}" "$L" "${SCHEME:-fe}" "${RTOL:-1.0e-4}" \
            "${EQUIL:-0.0}" "${PRE1D:-.false.}" .true. "$tol" -26000.0)"
    done
    ;;

  # The 3-D configuration at the truncations $LMAX_LIST_D. SUFFIX keeps the LGM
  # reference run (".ref") off the production ids.
  blockD-res)
    SCH=${SCHEME:-fe}; RT=${RTOL:-1.0e-4}
    for L in ${LMAX_LIST_D:-64 128 256}; do
      go "${MODEL_TAG:-D-V2}-3D.l$L.cold${SUFFIX:-}" \
         "$(nmlD "dF_3D_${L}${SUFFIX:-}" "$L" "$SCH" "$RT" "${EQUIL:-100000.0}" \
            "${PRE1D:-.true.}" .true. 1.0e-3 -26000.0 "${T1_D:-0.0}")"
    done
    ;;
esac

echo
echo "done. squeue -u $USER"
