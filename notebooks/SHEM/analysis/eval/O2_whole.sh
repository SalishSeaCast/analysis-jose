#!/bin/bash

#SBATCH --job-name=eval_O2_DFO
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=15:00:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-allen
# stdout and stderr file paths/names
#SBATCH --output=/scratch/jvalenti/run_SHEM/eval/stdout
#SBATCH --error=/scratch/jvalenti/run_SHEM/eval/stderr


RUN_ID="SHEM_eval_O2"
WORK_DIR="/home/jvalenti/MOAD/analysis-jose/notebooks/SHEM/analysis/eval"


cd ${WORK_DIR}
echo "working dir: $(pwd)"

echo "Starting run at $(date)"
pixi run -m $HOME/MOAD/analysis-jose -e shem python O2_whole.py
echo "Ended run at $(date)"
echo "Finished at $(date)"
