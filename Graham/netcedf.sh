#!/bin/bash

#SBATCH --job-name=Parcels_save
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=0
#SBATCH --time=5:59:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-allen
# stdout and stderr file paths/names
#SBATCH --output=/scratch/jvalenti/OParcels_runs/Parcels_outfall/results/stdout
#SBATCH --error=/scratch/jvalenti/OParcels_runs/Parcels_outfall/results/stderr


WORK_DIR="/scratch/jvalenti/OParcels_runs"

cd ${WORK_DIR}
echo "working dir: $(pwd)"

source /home/jvalenti/miniforge3/etc/profile.d/conda.sh
conda activate /home/jvalenti/miniforge3/envs/Parcels_v3

echo "Starting run at $(date)"
python3 -m zarr_tonetcdf.py
echo "Ended run at $(date)"
