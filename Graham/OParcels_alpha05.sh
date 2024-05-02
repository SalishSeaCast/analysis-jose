#!/bin/bash

#SBATCH --job-name=Parcel_alpha05
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=0
#SBATCH --time=47:59:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-allen
# stdout and stderr file paths/names
#SBATCH --output=/scratch/jvalenti/OParcels_runs/Parcels_alpha/results/stda05out
#SBATCH --error=/scratch/jvalenti/OParcels_runs/Parcels_alpha/results/stda05err


WORK_DIR="/scratch/jvalenti/OParcels_runs"

cd ${WORK_DIR}
echo "working dir: $(pwd)"

source /home/jvalenti/miniforge3/etc/profile.d/conda.sh
conda activate /home/jvalenti/miniforge3/envs/Parcels_24

echo "Starting run at $(date)"
mpirun -np 32 python3 -m Plastics_MCOP Parcels_alpha/Alpha05.yml
echo "Ended run at $(date)"
