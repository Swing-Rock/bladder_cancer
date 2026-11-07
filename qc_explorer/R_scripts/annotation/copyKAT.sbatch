#!/bin/bash
#SBATCH --job-name=copyKAT_qc_bladder_cancer
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 24:00:00
#SBATCH --output=%x_%j.out             # unique per job, written to the dir you sbatch from
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=YOUR_UCSD_EMAIL@ucsd.edu
#SBATCH --mail-type=all
#SBATCH -A CSD670


set -euo pipefail

# =============================================================================
# Edit the following before submission:
#
# 1. --mail-user
# 2. WORKDIR
# 3. RSCRIPT
# 4. SCRIPT
# 5. INPUT
# 6. SAMPLES (optional)
#
# =============================================================================

WORKDIR=/expanse/lustre/projects/csd670/<you>/qc_explorer  # the qc_explorer folder on Expanse
RSCRIPT=/home/<you>/miniconda3/envs/<env>/bin/Rscript      # Rscript from the env with copykat installed
SCRIPT=R_scripts/annotation/run_copykat.R                  # relative to WORKDIR or absolute
INPUT=datasets/<ID>/qc_output                              # QC output folder (or its samplesheet.csv)
SAMPLES=""                                                 # optional: only these samples, space-separated (blank = all)


cd "$WORKDIR"

# CPU thread settings
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
echo "Input:     $INPUT"
echo "Samples:   ${SAMPLES:-all}"
echo "CPUs:      ${SLURM_CPUS_PER_TASK:-NA}"
echo "============================================================"

"$RSCRIPT" "$SCRIPT" "$INPUT" $SAMPLES

echo "============================================================"
echo "End:       $(date)"
echo "JOB COMPLETE"
echo "============================================================"