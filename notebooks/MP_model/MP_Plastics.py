import sys 
import xarray as xr
import numpy as np
import os
import pandas as pd
import yaml
from numpy import random
import math
from datetime import datetime, timedelta
from parcels import FieldSet, Field, VectorField, ParticleSet, JITParticle, ParcelsRandom, Variable
from glob import glob
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter("ignore")

from Source.OP_functions import pandas_deploy

sys.path.append(os.getcwd()+'/Source')

from OP_Kernels import *
from OP_functions import *

def load_config(config_yaml):
   with open(config_yaml) as f:
       config = yaml.safe_load(f)
   return config

def MP_Plastics():
    param = load_config('yaml/config.yaml')
    start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
    length = param['param']['length'] # Set Time length [days] 
    dt = param['param']['dt'] #toggle between - or + to pick backwards or forwards 
    N = param['param']['N'] # Name of deploy loc file
    dd = param['param']['dd'] #max depth difference z in N
    name = param['file']['name'] #name output file
    dtp = param['param']['dtp'] #how often particle released in hours
    odt = param['param']['odt'] #how often data is recorded
    rrr = param['param']['r'] #radious of particle deployment
    SI = param['param']['SI'] #how many particles per super-individual

    duration = timedelta(days=length)
    print(f"The model will run for {duration.days} days, starting at {start}")

    lon, lat, z = pandas_deploy(N,SI,rrr,dd,dtp)
    N = len(lat)
    print(f"The model will release {N} particle every timestep")

    daterange = [start+timedelta(days=i) for i in range(length)]
    fn =  name + '.zarr'
    outfile = os.path.join(os.getcwd()+'/results/', fn)

    print(f"Output file: {outfile}")

    chunksize_fields = {
        "lon": ("lon", 66),   # matches 1 MPI box in lon
        "lat": ("lat", 28),   # matches 1 MPI box in lat
        "time": ("time", 1)
    }

    ####BUILD FIELDS FOR SIMULATION######
    #Fill in the list of variables that you want to use as fields
    varlist=['U','V','W']
    filenames,variables=filename_set(start,length,varlist)
    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw','time': 'time_counter'}
    field_set=FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=True)

    #Find file names and variable names ###'Diat','Flag'###
    varlist=['US','VS','WL','R','T','S','ssh','Bathy','Kz','totdepth']
    filenames,variables=filename_set(start,length,varlist)

    # #Add Stokes Drift fields
    # dimensions = {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
    # us = Field.from_netcdf(filenames['US'], variables['US'], dimensions,allow_time_extrapolation=True)
    # vs = Field.from_netcdf(filenames['VS'], variables['VS'], dimensions,allow_time_extrapolation=True)
    # wl = Field.from_netcdf(filenames['WL'], variables['WL'], dimensions,allow_time_extrapolation=True)
    # field_set.add_field(us)
    # field_set.add_field(vs)
    # field_set.add_field(wl)
    # field_set.add_vector_field(VectorField("stokes", us, vs, wl))

    #Add Vertical diffusivity coefficient field
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'depthw','time': 'time_counter'}
    Kz = Field.from_netcdf(filenames['Kz'], variables['Kz'], dimensions,allow_time_extrapolation=True)
    field_set.add_field(Kz)

    #Add fields located at node T
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'deptht','time': 'time_counter'}

    R = Field.from_netcdf(filenames['R'], variables['R'], dimensions,allow_time_extrapolation=True)
    S = Field.from_netcdf(filenames['S'], variables['S'], dimensions,allow_time_extrapolation=True)
    T = Field.from_netcdf(filenames['T'], variables['T'], dimensions,allow_time_extrapolation=True)
    field_set.add_field(R)
    field_set.add_field(S)
    field_set.add_field(T)

    #Add Bathymetry 2D field
    dimensions = {'lon': 'glamt', 'lat': 'gphit'}
    Bth = Field.from_netcdf(filenames['Bathy'], variables['Bathy'], dimensions,allow_time_extrapolation=True)
    TD = Field.from_netcdf(filenames['totdepth'], variables['totdepth'], dimensions,allow_time_extrapolation=True)
    field_set.add_field(Bth)
    field_set.add_field(TD)

    #Add SSH 
    dimensions = {'lon': 'glamt', 'lat': 'gphit','time': 'time_counter'}
    SSH = Field.from_netcdf(filenames['ssh'], variables['ssh'], dimensions,allow_time_extrapolation=True)
    field_set.add_field(SSH)

    #MPI decomposition
    Px = 6  # number of boxes along x (longitude)
    Py = 32   # number of boxes along y (latitude)\
    
    # Set domain bounds from the model grid
    xmin = field_set.U.grid.lon.min()
    xmax = field_set.U.grid.lon.max()
    ymin = field_set.U.grid.lat.min()
    ymax = field_set.U.grid.lat.max()

    x_edges = np.linspace(xmin, xmax, Px + 1)
    y_edges = np.linspace(ymin, ymax, Py + 1)
    def simple_partition_function(coords, mpi_size=1):
        """A very simple partition function
        that assigns particles to processors
        """
        return np.linspace(0, mpi_size, coords.shape[0], endpoint=False, dtype=np.int32)
    def nemo_partition_function(pset): 
        """Custom partition function for Parcels similar to NEMO MPI decomposition.
        Returns an array of processor IDs for each particle."""
        
        lon = pset.lon
        lat = pset.lat

        # Find which box each particle belongs to
        ix = np.searchsorted(x_edges, lon) - 1
        iy = np.searchsorted(y_edges, lat) - 1

        # Clamp indices to valid range
        ix = np.clip(ix, 0, Px - 1)
        iy = np.clip(iy, 0, Py - 1)

        # Map 2D box index to processor ID
        proc_id = iy * Px + ix

        return proc_id

    #Define our Parcels Particles
    class MPParticle(JITParticle):    
        diameter = Variable('diameter', initial = param['particle']['diameter'])
        length = Variable('length', initial = param['particle']['length'])
        Ub = Variable('Ub', initial = param['particle']['Ub'])  
        status = Variable('status', initial = 0)
        vvl_factor = Variable('fact', initial =  1)    
        ws = Variable('ws', initial =  0) 
        wa = Variable('wa', initial =  0) 
        wm = Variable('wm', initial =  0) 
        alpha = Variable('alpha',initial=param['particle']['alpha'])

    pset = ParticleSet.from_list(field_set, MPParticle, lon=lon, lat=lat, depth=z, repeatdt = timedelta(hours=dtp),partition_function=simple_partition_function)

    pset.execute([Advection,Buoyancy,turb_mix,Displacement,export,CheckOutOfBounds,KeepInOcean],
            runtime=duration, 
            dt=dt,
            output_file=pset.ParticleFile(name=outfile, outputdt=timedelta(hours=odt),chunks=(int(1e4), 1)))
    
if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print('Something went wrong')
    MP_Plastics()