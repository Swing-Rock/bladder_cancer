#!/bin/bash
#SBATCH --job-name=copykat_bladder
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G                      # shared partition charges max(cores, mem/~2GB) -- trim if you can
#SBATCH -t 24:00:00
#SBATCH --output=%x_%j.out             # unique per job, written to the dir you sbatch from
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=you@ucsd.edu       # <-- your real email, or delete this + the next line
#SBATCH --mail-type=END,FAIL
#SBATCH -A csd670

set -euo pipefail

# =============================================================================
# copyKAT on Expanse
# Fill in these three before submitting:
# =============================================================================
WORKDIR=/expanse/lustre/projects/csd670/<you>/uc_project   # project dir
RSCRIPT=/home/<you>/miniconda3/envs/<env>/bin/Rscript      # Rscript from the env with copykat installed
SCRIPT=scripts/run_copykat.R                               # relative to WORKDIR or absolute

# --- fail fast on unfilled / wrong paths ---
for v in WORKDIR RSCRIPT SCRIPT; do
  [[ -n "${!v}" ]] || { echo "ERROR: $v is empty" >&2; exit 1; }
done
cd "$WORKDIR"
[[ -x "$RSCRIPT" ]] || { echo "ERROR: $RSCRIPT not executable" >&2; exit 1; }
[[ -f "$SCRIPT"  ]] || { echo "ERROR: $SCRIPT not found in $WORKDIR" >&2; exit 1; }

# Threads: keep BLAS from oversubscribing on top of copyKAT's own parallelism
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
# Read this in R:  n.cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
export SLURM_CPUS_PER_TASK

echo "============================================================"
echo "copyKAT malignant cell labeling"
echo "Job ID:    ${SLURM_JOB_ID:-NA}"
echo "Node:      $(hostname)"
echo "Start:     $(date)"
echo "Workdir:   $(pwd)"
echo "Rscript:   $RSCRIPT"
echo "Script:    $SCRIPT"
echo "CPUs:      ${SLURM_CPUS_PER_TASK:-NA}"
echo "============================================================"

"$RSCRIPT" "$SCRIPT"

echo "============================================================"
echo "End:       $(date)"
echo "JOB COMPLETE"
echo "============================================================"
