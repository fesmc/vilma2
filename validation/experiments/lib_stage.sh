# lib_stage.sh — machine selection and job staging. Source, do not execute.
#
# Defines VALIDATION (this directory's parent), MODEL (the repository root) and
# the V_* settings of machines/<name>.sh. The machine is VALIDATION_MACHINE, or
# is taken from the host name.

VALIDATION="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
MODEL="$(cd "$VALIDATION/.." && pwd -P)"
export PATH="$HOME/.local/bin:$PATH"          # runme

v_machine() {
  if [ -n "${VALIDATION_MACHINE:-}" ]; then echo "$VALIDATION_MACHINE"; return; fi
  case "$(hostname -f 2>/dev/null || hostname)" in
    albedo*|prod-*|fat-*|gpu-*) echo albedo ;;
    levante*|l[0-9]*)           echo levante ;;
    *) echo "lib_stage: unknown host $(hostname); set VALIDATION_MACHINE to one of: $(cd "$VALIDATION/machines" && ls | sed 's/\.sh$//' | tr '\n' ' ')" >&2
       return 1 ;;
  esac
}
V_MACHINE="$(v_machine)" || exit 1
[ -f "$VALIDATION/machines/$V_MACHINE.sh" ] || { echo "lib_stage: no machines/$V_MACHINE.sh" >&2; exit 1; }
# shellcheck disable=SC1090
. "$VALIDATION/machines/$V_MACHINE.sh"

# v_busy <dir>: true when one of the user's queued or running jobs works in <dir>.
# runme clears a directory it stages into, and a running job then dies at its
# next write to out.nc.
v_busy() { squeue -u "$USER" -h -o "%Z" 2>/dev/null | grep -xF "$1" >/dev/null; }

# v_fix_job <job.submit> <threads> [bind]
# Rewrites runme's job script so the thread count is set once, before srun.
# V_EXCLUSIVE=1: whole node and all its memory. V_EXCLUSIVE=0: as many cores as
# threads. After srun the job prints its id, node and elapsed time to out.out,
# in the form analysis/BlockA.jl reads (slurm_meta).
v_fix_job() {
  python3 - "$1" "$2" "${3:-spread}" "$V_NCORE" "$V_EXCLUSIVE" <<'PY'
import re, sys
path, omp, bind, ncore, excl = sys.argv[1:6]
s = open(path).read()
cpus = ncore if excl == '1' else omp
if excl == '1' and '--exclusive' not in s:
    s = s.replace('#SBATCH --partition=',
                  '#SBATCH --exclusive\n#SBATCH --mem=0\n#SBATCH --partition=', 1)
s = re.sub(r'^#SBATCH --cpus-per-task=\d+', f'#SBATCH --cpus-per-task={cpus}', s, flags=re.M)
# runme writes the thread count more than once; keep one block before srun.
s = re.sub(r'^export (OMP_NUM_THREADS|SRUN_CPUS_PER_TASK|MKL_NUM_THREADS|OMP_PROC_BIND|OMP_PLACES)=.*\n',
           '', s, flags=re.M)
block = (f'export OMP_NUM_THREADS={omp}\nexport MKL_NUM_THREADS={omp}\n'
         f'export SRUN_CPUS_PER_TASK={cpus}\nexport OMP_PROC_BIND={bind}\n'
         'export OMP_PLACES=cores\n')
tail = (r"""
rc=$?
printf '* JobID            : %s\n' "$SLURM_JOB_ID"
printf '* Nodelist         : %s\n' "$SLURM_JOB_NODELIST"
printf '* Elapsed time     : %02d:%02d:%02d\n' $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60))
exit $rc
""")
s, n = re.subn(r'^(srun .*)$', lambda m: block + m.group(1) + tail, s, count=1, flags=re.M)
if n != 1:
    sys.exit(f'v_fix_job: no srun line in {path}')
open(path, 'w').write(s)
PY
}

# v_stage <dir> <exe> <nml> <wall> <threads> [bind]
# Stages a run directory with runme (binary, namelist, data and input links) and
# fixes its job script. Does not submit.
v_stage() {
  local dir=$1 exe=$2 par=$3 wall=$4 omp=$5 bind=${6:-spread}
  if v_busy "$dir"; then
    echo "REFUSING: a job is queued or running in $dir (cancel it, or stage under another id)" >&2
    return 1
  fi
  ( cd "$MODEL" && runme -o "$dir" -e "$exe" -n "$par" \
        -q "$(v_queue "$wall")" -w "$wall" --omp "$omp" -s ) >/dev/null
  [ -f "$dir/job.submit" ] || { echo "v_stage: runme did not stage $dir (try: cd $MODEL && runme info)" >&2; return 1; }
  v_fix_job "$dir/job.submit" "$omp" "$bind"
}
