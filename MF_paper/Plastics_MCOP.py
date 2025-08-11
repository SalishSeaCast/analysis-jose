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


sys.path.append('/home/jvalenti/MOAD/analysis-jose/Graham/Source')
from OP_Kernels import *

def load_config(config_yaml):
   with open(config_yaml[0]) as f:
       config = yaml.safe_load(f)
   return config

def zarr_tonet(fileoutname):
    from os import path
    name =  fileoutname.split('_')[1]
    fname = fileoutname
    print(fileoutname)
    files = glob(path.join(fname, "proc*"))
    ds = xr.concat(
        [xr.open_zarr(f) for f in files],
        dim="trajectory",
        compat="no_conflicts",
        coords="minimal",
    )
    return ds,name

def simple_partition_function(coords, mpi_size=1):
    """A very simple partition function
    that assigns particles to processors
    """
    return np.linspace(0, mpi_size, coords.shape[0], endpoint=False, dtype=np.int32)

def pandas_deploy(N,MFc,r,dd,dtp):
    n =1
    MFc = float(MFc)
    Outfall_deploy = pd.read_csv(N, index_col = [0])
    Pol = list(Outfall_deploy.Population)
    Lat = Outfall_deploy.Lat
    Lon = Outfall_deploy.Lon
    Depth = Outfall_deploy.Depth
    clat = []
    clon = []
    cz = []
    for i,loc in enumerate(Pol):
        for j in range(int(round((loc*250*dtp)/MFc,0))):
            clat.append(Lat.iat[i])
            clon.append(Lon.iat[i])
            try:
                cz.append(float(Depth.iat[i]))
            except ValueError:
                cz.append(2)
    N = len(clat)
    deg2m = 111000 * np.cos(49 * np.pi / 180)
    var = (r / (deg2m * 3))**2
    x_offset, y_offset = np.random.multivariate_normal([0, 0], [[var, 0], [0, var]], [n,N]).T
    z_offset=np.random.random_sample([N]).T*(dd)
    lon = np.zeros([N])
    lat = np.zeros([N])
    z = np.zeros([N])
    for i in range(N):
        lon[i]=(clon[i] + x_offset[i])
        lat[i]=(clat[i] + y_offset[i])
        z[i]=(cz[i] + z_offset[i])
    return lon,lat,z

def filename_set(start,length,varlist=['U','V','W']):
    '''filename,variables,dimensions = filename_set(start,duration,varlist=['U','V','W'],local=1)
    Modify function to include more default variables
    define start as: e.g, datetime(2018, 1, 17)
    length= number of days'''
    
    duration = timedelta(days=length)
    #Build filenames
    paths = path()
    Tlist,Ulist, Vlist, Wlist = [], [], [], [], []
    Waveslist = []
   
    for day in range(duration.days):
        path_NEMO = make_prefix(start + timedelta(days=day), paths['NEMO'])
        path_NEMO_d = make_prefix(start + timedelta(days=day), paths['NEMO'],res='d')
        Ulist.append(path_NEMO + '_grid_U.nc')
        Vlist.append(path_NEMO + '_grid_V.nc')
        Wlist.append(path_NEMO + '_grid_W.nc')
        Tlist.append(path_NEMO + '_grid_T.nc')
        Waveslist.append(get_WW3_path(start + timedelta(days=day)))        

    # Load NEMO forcing 
    filenames = {
        'U': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Ulist},
        'V': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Vlist},
        'W': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Wlist},
        'Kz': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Wlist},
        'T': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Tlist[0], 'data': Tlist},
        'S': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Tlist[0], 'data': Tlist},
        'ssh': {'lon': paths['coords'], 'lat': paths['coords'], 'data': Tlist},
        'R': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Tlist[0], 'data': Tlist},
        'Bathy' : {'lon': paths['coords'], 'lat': paths['coords'], 'data': paths['bat']},
        'gdepth' : {'lon': paths['coords'], 'lat': paths['coords'],'depth': Wlist[0], 'data': paths['mask']},
        'totdepth' : {'lon': paths['coords'], 'lat': paths['coords'], 'data': paths['mask']},
        'US' : {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist},
        'VS' : {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist},
        'WL' : {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist}
    }
    
    variables = {'U': 'vozocrtx', 'V': 'vomecrty','W': 'vovecrtz','T':'votemper','S':'vosaline','R':'sigma_theta',
        'US':'uuss','VS':'vuss','WL':'lm','Bathy':'Bathymetry','Kz':'vert_eddy_diff','ssh':'sossheig','totdepth':'totaldepth'}
        
    file2,var2 = {},{}
    for var in varlist:
        file2[var]=filenames[var]
        var2[var]=variables[var]
    return file2,var2

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
    MFc = param['param']['MFc']

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


    ds,name = zarr_tonet(outfile)
    ds.to_netcdf('/home/jvalenti/projects/rrg-allen/jvalenti/'+name+'.nc')

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print('Something went wrong')
    Plastics_OP(config)
     
