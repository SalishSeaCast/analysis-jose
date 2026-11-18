# bench.py <latlon> <depth>
import sys
from datetime import timedelta
from mpi_Parcels import mpi_Parcels

ll, dd = int(sys.argv[1]), int(sys.argv[2])
cs = {"time": ("time", 1), "depth": ("depth", 40),
      "lat": ("latitude", 200), "lon": ("longitude", 400)}   # same convention as the yaml

mpi_Parcels('yaml/config.yaml', runtime=timedelta(hours=6), outputdt=timedelta(hours=1),
            outfile=f'/home/jvalenti/scratch/Parcels/results/bench_{ll}.zarr', cs_override='none', convert=False, seed=1,overwrite=True)