#!/bin/bash

#SBATCH --job-name=eval_O2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=8000
#SBATCH --time=4:00:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-allen
# stdout and stderr file paths/names
#SBATCH --output=/scratch/jvalenti/run_SHEM/eval/stdout
#SBATCH --error=/scratch/jvalenti/run_SHEM/eval/stderr


RUN_ID="SHEM_eval_O2_mort_hbac"
WORK_DIR="/home/jvalenti/MOAD/analysis-jose/notebooks/SHEM/tuning"


cd ${WORK_DIR}
echo "working dir: $(pwd)"

export PATH=/home/jvalenti/miniforge3/envs/SHEM/bin:$PATH
source /home/jvalenti/miniforge3/etc/profile.d/conda.sh
conda activate SHEM

echo "Starting run at $(date)"
python O2_puget.py SSBase
echo "Ended run at $(date)"
echo "Finished at $(date)"
