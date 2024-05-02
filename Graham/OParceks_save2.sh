#!/bin/bash

#SBATCH --job-name=Saving_tonet3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=0
#SBATCH --time=05:59:00
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL
#SBATCH --account=rrg-allen
# stdout and stderr file paths/names
#SBATCH --output=/home/jvalenti/MOAD/analysis-jose/Graham/results/stdsaveout
#SBATCH --error=/home/jvalenti/MOAD/analysis-jose/Graham/results/stdsaveerr


WORK_DIR="/home/jvalenti/MOAD/analysis-jose/Graham"

cd ${WORK_DIR}
echo "working dir: $(pwd)"

source /home/jvalenti/miniforge3/etc/profile.d/conda.sh
conda activate /home/jvalenti/miniforge3/envs/Parcels_24

echo "Starting save at $(date)"
python3 -m zarr_tonetcdf Run_Julio20190701_1n_20190729_1n.zarr
echo "Ended save at $(date)"
