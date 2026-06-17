#!/bin/bash

#SBATCH --job-name=eval_NO3_PS
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

CONFIGS=(SSBase SHEM18 diat_pref exc_hbac exc_hbac_2 growth_flag growth_flag_2 mort_hbac mort_hbac_2 pred_flag remin remin2 remin2_l predmine predmine_z2)

CONFIG=${CONFIGS[$SLURM_ARRAY_TASK_ID]}

RUN_ID="SHEM_eval_NO3_puget"
WORK_DIR="/home/jvalenti/MOAD/analysis-jose/notebooks/SHEM/tuning"

cd ${WORK_DIR}
echo "working dir: $(pwd)"

echo "Running config: ${CONFIG}"
echo "Starting run at $(date)"

pixi run -e shem python NO3_puget_v2.py ${CONFIG}

echo "Ended run at $(date)"
