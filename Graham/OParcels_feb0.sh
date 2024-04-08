#!/bin/bash

#SBATCH --job-name=Parcel_Febrero
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=0
#SBATCH --time=23:59:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-allen
# stdout and stderr file paths/names
#SBATCH --output=/scratch/jvalenti/OParcels_runs/Parcels_alpha/results/stdafebout
#SBATCH --error=/scratch/jvalenti/OParcels_runs/Parcels_alpha/results/stdafeberr


WORK_DIR="/home/jvalenti/MOAD/analysis-jose/Graham"

cd ${WORK_DIR}
echo "working dir: $(pwd)"

source /home/jvalenti/miniforge3/etc/profile.d/conda.sh
conda activate /home/jvalenti/miniforge3/envs/Parcels_24

echo "Starting run at $(date)"
mpirun -np 32 python3 -m Plastics_MCOP yaml/Feb0.yml
echo "Ended run at $(date)"
