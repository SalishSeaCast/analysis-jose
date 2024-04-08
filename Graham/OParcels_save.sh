#!/bin/bash

#SBATCH --job-name=Saving_tonet1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#BATCH --mem=32G
#SBATCH --time=05:59:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-allen
# stdout and stderr file paths/names
#SBATCH --output=/scratch/jvalenti/OParcels_runs/Parcels_outfall/results/stdsaveout
#SBATCH --error=/scratch/jvalenti/OParcels_runs/Parcels_outfall/results/stdsaveerr


WORK_DIR="/scratch/jvalenti/OParcels_runs"

cd ${WORK_DIR}
echo "working dir: $(pwd)"

source /home/jvalenti/miniforge3/etc/profile.d/conda.sh
conda activate /home/jvalenti/miniforge3/envs/Parcels_24

echo "Starting save at $(date)"
python3 -m zarr_tonetcdf Runs_alpha120190101_1n_20190201_1n.zarr
echo "Ended save at $(date)"
