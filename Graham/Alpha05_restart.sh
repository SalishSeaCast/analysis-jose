#!/bin/bash

#SBATCH --job-name=Restart_alpharun05
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=16G
#SBATCH --time=23:59:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-allen
# stdout and stderr file paths/names
#SBATCH --output=/home/jvalenti/MOAD/analysis-jose/Graham/results/stdR05out
#SBATCH --error=/home/jvalenti/MOAD/analysis-jose/Graham/results/stdR05err


WORK_DIR="/scratch/OParcels_runs/MOAD/analysis-jose/Graham" 

cd ${WORK_DIR}
echo "working dir: $(pwd)"

source /home/jvalenti/miniforge3/etc/profile.d/conda.sh
conda activate /home/jvalenti/miniforge3/envs/Parcels_24

echo "Starting run at $(date)"
mpirun -np 32 python3 -m Plastics_MCrestart yaml/alpha05_restart.yaml
echo "Ended run at $(date)"
