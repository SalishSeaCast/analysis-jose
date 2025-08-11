import sys 
import xarray as xr
import numpy as np
import os
import pandas as pd
import yaml
from numpy import random
#from random import randint
import math
from datetime import datetime, timedelta
from parcels import FieldSet, Field, VectorField, ParticleSet, JITParticle, ParcelsRandom, Variable
from glob import glob


sys.path.append('/Source')

from OP_Kernels import *
from OP_functions import *


def Plastics_OP(config):

    param = load_config(config)

    start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
    length = param['param']['length'] # Set Time length [days] 
    dt = param['param']['dt'] #toggle between - or + to pick backwards or forwards 
    N = param['param']['N'] # Name of deploy loc file
    dd = param['param']['dd'] #max depth difference z in N
    name = param['file']['name'] #name output file
    dtp = param['param']['dtp'] #how often particle released in hours
    odt = param['param']['odt'] #how often data is recorded
    rrr = param['param']['r'] #radious of particle deployment
    MFc = param['param']['MFc'] #how many particles per super-individual

    duration = timedelta(days=length)
    print(f"The model will run for {duration.days} days, starting at {start}")

    lon, lat, z = pandas_deploy(N,MFc,rrr,dd,dtp)
    N = len(lat)
    print(f"The model will release {N} particle every timestep")

    daterange = [start+timedelta(days=i) for i in range(length)]
    fn =  name + '.zarr'
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
    
    R = Field.from_netcdf(filenames['R'], variables['R'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    S = Field.from_netcdf(filenames['S'], variables['S'], dimensions,allow_time_extrapolation=True, chunksize='auto')
    T = Field.from_netcdf(filenames['T'], variables['T'], dimensions,allow_time_extrapolation=True, chunksize='auto')
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

    #Define our Ocean Parcels Particles
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

    ######RUN OCEAN PARCELS WITH DEFINED PARTICLE AND PRESET FIELDS
    if dtp == 0:
        pset = ParticleSet.from_list(field_set, MPParticle, lon=lon, lat=lat, depth=z,time=start+timedelta(hours=odt),partition_function=simple_partition_function)
    else:
        pset = ParticleSet.from_list(field_set, MPParticle, lon=lon, lat=lat, depth=z, repeatdt = timedelta(hours=dtp),partition_function=simple_partition_function)

    pset.execute([Advection,Buoyancy,Stokes_drift,turb_mix,Displacement,export,CheckOutOfBounds,KeepInOcean],
        runtime=duration, 
        dt=dt,
        output_file=pset.ParticleFile(name=outfile, outputdt=timedelta(hours=odt),chunks=(int(1e4), 1)))

    #Convert and save output as netcdf 
    ds,name = zarr_tonet(outfile)
    ds.to_netcdf('/home/jvalenti/projects/rrg-allen/jvalenti/'+name+'.nc')

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print('Something went wrong')
    Plastics_OP(config)
     
