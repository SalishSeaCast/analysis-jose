#!/bin/bash
#SBATCH --job-name=shem_2023
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=0
#SBATCH --time=6:00:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=def-allen
# stdout and stderr file paths/names
#SBATCH --output=/scratch/jvalenti/stdout
#SBATCH --error=/scratch/jvalenti/stderr

source ~/.bashrc

WORK_DIR="/home/jvalenti/MOAD/analysis-jose/notebooks/SHEM/analysis"
cd ${WORK_DIR}


pixi run -e shem python dayproc.py