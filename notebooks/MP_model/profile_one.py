# profile_one.py
import sys, cProfile, pstats
import numpy as np
from datetime import timedelta
sys.path.append('/home/jvalenti/MOAD/analysis-jose/notebooks/MP_model/Source')

from OP_functions import load_config, particle_source
from mpi_Parcels import mpi_Parcels

cfg = 'yaml/config.yaml'
param = load_config(cfg)

MP = particle_source(param, seed=1234, runtime_hours=1)      # all released in the first hour
rng = np.random.default_rng(0)
idx = np.sort(rng.choice(len(MP.lat), 1000, replace=False))
particles = {k: np.asarray(v)[idx] for k, v in vars(MP).items()}

cs = {"time": ("time", 1), "depth": ("depth", 40),
      "lat": ("latitude", 200), "lon": ("longitude", 400)}

def go():
    mpi_Parcels(cfg, runtime=timedelta(hours=6), outputdt=timedelta(hours=1),
                outfile='/home/jvalenti/scratch/Parcels/bench/prof_one.zarr',
                particles=particles, cs_override=cs, convert=False,
                seed=1, overwrite=True)

cProfile.run('go()', 'prof_one.stats')
p = pstats.Stats('prof_one.stats')
p.sort_stats('cumulative').print_stats(40)