#!/bin/bash
#SBATCH --job-name=copyKAT_qc
#SBATCH --output=output.txt
#SBATCH -p shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --export=ALL
#SBATCH -t 24:00:00
#SBATCH --mail-user=YOUR_EMAIL_HERE
#SBATCH --mail-type=all
#SBATCH -A CSD670
#SBATCH --mem=100G

module load singularitypro

echo "Job started"

SIF=/expanse/lustre/projects/csd670/akao1/cell2location/cell2loc.sif
DATA_DIR=/PATH/TO/PROJECT/DIRECTORY/THAT/YOU/MADE
SCRIPT=$DATA_DIR/YOUR_CELL2LOCATION_SCRIPT.py

singularity exec --bind $DATA_DIR:$DATA_DIR $SIF which python
singularity exec --bind $DATA_DIR:$DATA_DIR $SIF python --version
singularity exec --bind $DATA_DIR:$DATA_DIR $SIF python $SCRIPT

echo "Job finished"