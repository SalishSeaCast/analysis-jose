# bench_stage.py <tag>
import sys, os
from datetime import timedelta
sys.path.append('/home/jvalenti/MOAD/analysis-jose/notebooks/MP_model/Source')

from OP_functions import load_config, particle_source
from mpi_Parcels import mpi_Parcels

tag = sys.argv[1]                      # "scratch" or "staged"
cfg = 'yaml/config.yaml'
param = load_config(cfg)
MP = particle_source(param, seed=1234, runtime_hours=12)    # N = 17,005, fixed

mpi_Parcels(cfg, runtime=timedelta(hours=12), outputdt=timedelta(hours=1),
            outfile=f'/home/jvalenti/scratch/Parcels/bench/stage_{tag}.zarr',
            particles=vars(MP), cs_override='none', convert=False, seed=1, overwrite=True)