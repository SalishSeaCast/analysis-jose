#!/bin/bash 
#SBATCH --job-name=parcels_bench
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --account=def-allen    
#SBATCH --mail-user=jvalenti@eoas.ubc.ca
#SBATCH --mail-type=ALL

# stdout and stderr file paths/names
#SBATCH --output=/home/jvalenti/scratch/Parcels/results/stdout
#SBATCH --error=/home/jvalenti/scratch/Parcels/results/stderr

module load StdEnv/2023
module load netcdf-fortran-mpi/4.6.1

export OMP_NUM_THREADS=1
export TMPDIR=$SLURM_TMPDIR          # compile kernels on node-local storage

cd /home/jvalenti/MOAD/analysis-jose/notebooks/MP_model   # where bench.py and yaml/ live

module unload openmpi/4.1.5
module unload ucx ucc libfabric pmix 2>/dev/null     # optional, silences UCX noise
module list 2>&1 | grep -q 'gcc/12.3' || module load gcc/12.3   # the JIT needs gcc

for v in $(env | grep -oE '^(PMIX|PMI|OMPI)_[A-Za-z0-9_]*'); do unset $v; done

for n in 35 45 55; do
    pixi run -e parcels mpirun -np $n python bench_stage.py 200 40 | grep BENCH
done