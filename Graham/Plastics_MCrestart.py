import sys 
import xarray as xr
import numpy as np
import os
import yaml
from numpy import random
#from random import randint
import math
import warnings
from datetime import datetime, timedelta
from parcels import FieldSet, Field, VectorField, ParticleSet, JITParticle, ParcelsRandom, Variable
warnings.simplefilter(action='ignore', category=FutureWarning)

sys.path.append('/home/jvalenti/MOAD/analysis-jose/Graham/Source')
from OP_functions import *

def simple_partition_function(coords, mpi_size=1):
    """A very simple partition function
    that assigns particles to processors
    """
    return np.linspace(0, mpi_size, coords.shape[0], endpoint=False, dtype=np.int32)
def Restart_OP(config):
    local = 0
    param = load_config(config)
#Definitions
    start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
    length = param['param']['length'] # Set Time length [days] 
    dt = param['param']['dt'] #toggle between - or + to pick backwards or forwards 
    restart = param['param']['N'] # Name of deploy loc file
    dd = param['param']['dd'] #max depth difference z in N
    name = param['file']['name'] #name output file
    dtp = param['param']['dtp'] #how often particle released in hours
    odt = param['param']['odt'] #how often data is recorded
    rrr = param['param']['r'] #radious of particle deployment
    MFc = param['param']['MFc']

    duration = timedelta(days=length)
    print(f"The model will run for {duration.days} days, starting at {start}")

    daterange = [start+timedelta(days=i) for i in range(length)]

    fn =  name +'_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.zarr'   
    outfile = os.path.join('/scratch/jvalenti/OParcels_runs/Parcels_alpha/results/', fn)

####BUILD FIELDS FOR SIMULATION######
    #Fill in the list of variables that you want to use as fields
    varlist=['U','V','W']
    filenames,variables=filename_set(start,length,varlist)
    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw','time': 'time_counter'}
    field_set=FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=True, chunksize='auto')

    #Find file names and variable names ###'Diat','Flag'###
    varlist=['US','VS','WL','R','T','S','ssh','Bathy','Kz','totdepth','Vol']
    filenames,variables=filename_set(start,length,varlist)

    #Add Stokes Drift fields
    dimensions = {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
    us = Field.from_netcdf(filenames['US'], variables['US'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    vs = Field.from_netcdf(filenames['VS'], variables['VS'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    wl = Field.from_netcdf(filenames['WL'], variables['WL'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(us)
    field_set.add_field(vs)
    field_set.add_field(wl)
    field_set.add_vector_field(VectorField("stokes", us, vs, wl))

    #Add Vertical diffusivity coefficient field
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'depthw','time': 'time_counter'}
    Kz = Field.from_netcdf(filenames['Kz'], variables['Kz'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(Kz)

    #Add fields located at node T
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'deptht','time': 'time_counter'}
    #Diat = Field.from_netcdf(filenames['Diat'], variables['Diat'], dimensions,allow_time_extrapolation=True)
    #Flag = Field.from_netcdf(filenames['Flag'], variables['Flag'], dimensions,allow_time_extrapolation=True)
    R = Field.from_netcdf(filenames['R'], variables['R'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    S = Field.from_netcdf(filenames['S'], variables['S'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    T = Field.from_netcdf(filenames['T'], variables['T'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    #field_set.add_field(Diat)
    #field_set.add_field(Flag)
    field_set.add_field(R)
    field_set.add_field(S)
    field_set.add_field(T)

    #Add Bathymetry 2D field
    dimensions = {'lon': 'glamt', 'lat': 'gphit'}
    Bth = Field.from_netcdf(filenames['Bathy'], variables['Bathy'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    TD = Field.from_netcdf(filenames['totdepth'], variables['totdepth'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(Bth)
    field_set.add_field(TD)

    #Add Volume 3D field
    dimensions = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'depthw'}
    Vol = Field.from_netcdf(filenames['Vol'], variables['Vol'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(Vol)

    #Add SSH 
    dimensions = {'lon': 'glamt', 'lat': 'gphit','time': 'time_counter'}
    SSH = Field.from_netcdf(filenames['ssh'], variables['ssh'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    field_set.add_field(SSH)
    #Fraser = Field.from_netcdf(filenames['FS'], variables['FS'], dimensions,allow_time_extrapolation=True,timestamps=get_timestamps(start,length))
    #field_set.add_field(Fraser) 
    class MPParticle(JITParticle):    
        diameter = Variable('diameter', initial = param['particle']['diameter'])
        length = Variable('length', initial = param['particle']['length'])
        Ub = Variable('Ub', initial = param['particle']['Ub'])  #days to have 67% probability of unbeaching
        tau = Variable('tau', initial =  param['particle']['tau']) # track number of particles
        status = Variable('status', initial = 0)
        fratio = Variable('fratio', initial = param['particle']['fratio'])
        fact = Variable('fact', initial =  1)    
        ws = Variable('ws', initial =  0) 
        wa = Variable('wa', initial =  0) 
        wm = Variable('wm', initial =  0) 
        cellvol = Variable('cellvol', initial =  0) # MF per parcel
        alpha = Variable('alpha',initial=param['particle']['alpha'])

    ds = xr.open_dataset(restart,decode_times=False)
    st = int(ds.time[0,-1].values)
    year,month,day = ds.time.units.split(' ')[2].split('-')
    day = day.split('T')[0]
    restartT = datetime(int(year), int(month), int(day))+ timedelta(seconds= st)
    DS = ds.to_dataframe()
    dslast = DS[DS.time==st]
    dslast = dslast[dslast.status == 1]
    
######RUN OCEAN PARCELS WITH DEFINED PARTICLE AND PRESET FIELDS 
    pset = ParticleSet.from_list(field_set, MPParticle, lon=dslast.lon, lat=dslast.lat, depth=dslast.z,time=restartT+timedelta(hours=odt),partition_function=simple_partition_function)
    print('load particles successful')
    pset.execute([Advection,Buoyancy,Stokes_drift,turb_mix,Displacement,export,CheckOutOfBounds,KeepInOcean],
        runtime=duration, 
        dt=dt,
        output_file=pset.ParticleFile(name=outfile, outputdt=timedelta(hours=odt),chunks=(int(1e4), 1)))
    
    from glob import glob
    from os import path
    names =  outfile.split('_')[1]
    files = glob(path.join(outfile, "proc*"))
    ds = xr.concat(
        [xr.open_zarr(f) for f in files],
        dim="trajectory",
        compat="no_conflicts",
        coords="minimal",
    )
    ds.to_netcdf('/home/jvalenti/projects/rrg-allen/jvalenti/restart'+names+'.nc')


def get_timestamps(start,length):
    timestamps=[]
    duration = timedelta(days=length)
    for day in range(duration.days):
        timestamps.append([start + timedelta(days=day)])
    return np.array(timestamps, dtype='datetime64')

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except ValueError:
            print('Something went wrong')
    Restart_OP(config)
     
