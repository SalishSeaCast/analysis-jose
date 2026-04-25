#!/bin/bash

#SBATCH --job-name=eval_O2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=8000
#SBATCH --time=3:00:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-allen
# stdout and stderr file paths/names
#SBATCH --output=/scratch/jvalenti/run_SHEM/eval/stdout
#SBATCH --error=/scratch/jvalenti/run_SHEM/eval/stderr


RUN_ID="SHEM_eval_O2_predmine2l_z2"
WORK_DIR="/home/jvalenti/MOAD/analysis-jose/notebooks/SHEM/tuning"


cd ${WORK_DIR}
echo "working dir: $(pwd)"

echo "Starting run at $(date)"
pixi run -e shem python O2_puget_v2.py predmine2l_z2
echo "Ended run at $(date)"
echo "Finished at $(date)"
