import sys
import xarray as xr
import numpy as np
import os
import yaml
import math
from datetime import datetime, timedelta
from parcels import FieldSet, Field, VectorField, ParticleSet, JITParticle, ErrorCode, ParcelsRandom, Variable

sys.path.append('/home/jvalenti/MOAD/analysis-jose/notebooks/parcels')
from Kernels_biofilm import DeleteParticle, Buoyancy, AdvectionRK4_3D, Stokes_drift, Beaching, Unbeaching
from OP_functions_biofilm import *

def fibers_OP(config,local=0,restart=0):
    param = load_config(config)
    #Definitions
    start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
    Tmax = param['param']['length'] # Set Time length [days] 
    dt = param['param']['dt'] #toggle between - or + to pick backwards or forwards 
    N = param['param']['N'] # number of deploying locations
    n = param['param']['n'] # 1000   # number of particles per location
    dmin = param['param']['dmin'] #minimum depth
    dd = param['param']['dd'] #max depth difference from dmin
    name = param['file']['name'] #name output file
    dtp = param['param']['dtp'] #how often particle released in hours
# Define paths
    paths = path(local)
#Set outfall coordinates (Modify to choose other deploying location)    
    coord=xr.open_dataset(paths['coords'],decode_times=False)
    outf_lat=coord['nav_lat'][445,304]
    outf_lon=coord['nav_lon'][445,304]
    clon, clat = [float(outf_lon)],[float(outf_lat)] 

    duration = timedelta(days=Tmax)
    x_offset, y_offset, z = p_deploy(N,n,dmin,dd)

#Set deploy locations
    lon = np.zeros([N,n])
    lat = np.zeros([N,n])
    for i in range(N):
        lon[i,:]=(clon[i] + x_offset[i,:])
        lat[i,:]=(clat[i] + y_offset[i,:])

#Set start date time and the name of the output file

    daterange = [start+timedelta(days=i) for i in range(Tmax)]
    fn =  name + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.nc'
    outfile = os.path.join(paths['out'], fn)

####BUILD FIELDS FOR SIMULATION######

#Fill in the list of variables that you want to use as fields
    varlist=['U','V','W','R']
    filenames,variables,dimensions=filename_set(start,Tmax,varlist,local)
    field_set=FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=True)

    varlist=['US','VS','WL']
    filenames,variables,dimensions=filename_set(start,Tmax,varlist,local)

    us = Field.from_netcdf(filenames['US'], variables['US'], dimensions,allow_time_extrapolation=True)
    vs = Field.from_netcdf(filenames['VS'], variables['VS'], dimensions,allow_time_extrapolation=True)
    wl = Field.from_netcdf(filenames['WL'], variables['WL'], dimensions,allow_time_extrapolation=True)
    field_set.add_field(us)
    field_set.add_field(vs)
    field_set.add_field(wl)
    field_set.add_vector_field(VectorField("stokes", us, vs, wl))

    filenames,variables,dimensions=filename_set(start,Tmax,['Bathy'],local)
    Bth = Field.from_netcdf(filenames['Bathy'], variables['Bathy'], dimensions,allow_time_extrapolation=True)
    field_set.add_field(Bth)
    MPParticle = particle_maker(param)
    
######RUN OCEAN PARCELS WITH DEFINED PARTICLE AND PRESET FIELDS
    if restart==1:
        name_temp=find_temp(paths['out'])
        os.system(f"cd {paths['out']} && parcels_convert_npydir_to_netcdf {name_temp}")
        outfile=newest(paths['out'])
        
        pset = ParticleSet.from_particlefile(field_set, MPParticle,outfile)
    else:
        pset = ParticleSet.from_list(field_set, MPParticle, lon=lon, lat=lat, depth=z, repeatdt = timedelta(hours=dtp))
    
    k_sink = pset.Kernel(Buoyancy)
    k_waves = pset.Kernel(Stokes_drift)
    k_beach = pset.Kernel(Beaching)
    k_unbeach = pset.Kernel(Unbeaching)
    
    pset.execute(AdvectionRK4_3D + k_sink + k_waves + k_beach + k_unbeach,
                runtime=duration, 
                dt=dt,
                output_file=pset.ParticleFile(name=outfile, outputdt=timedelta(hours=1)),
                recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

def load_config(config_yaml):
   with open(config_yaml[0]) as f:
       config = yaml.safe_load(f)
   return config

def particle_maker(config):
    #Define particle properties 
    class MPParticle(JITParticle):
        if 'ro' in config['particle']:        
            ro = Variable('ro', initial = config['particle']['ro'])  # config['particle']['ro']
        if 'diameter' in config['particle']:           
            diameter = Variable('diameter', initial = config['particle']['diameter'])
        if 'length' in config['particle']:  
            length = Variable('length', initial = config['particle']['length'])
        if 'Lb' in config['particle']:  
            Lb = Variable('Lb', initial = config['particle']['Lb'])  #days needed in days for particle to have 67% probability of beaching if in beaching zone (500m)
        if 'Db' in config['particle']:  
            Db = Variable('Db', initial = config['particle']['Db']) #Distance at which particles can randomly beach.
        if 'Ub' in config['particle']:  
            Ub = Variable('Ub', initial = config['particle']['Ub'])  #days to have 67% probability of unbeaching
        if 'beached' in config['particle']:  
            beached = Variable('beached', initial = 0)
        if 'Ws' in config['particle']:  
            Ws = Variable('Ws', initial =  config['particle']['Ws']) #200m/dia
        if 'tau' in config['particle']:  
            tau = Variable('tau', initial =  0) # track age particle
    return MPParticle

def find_temp(rootdir):
    dirs=[]
    for file in os.listdir(rootdir):
        d = os.path.join(rootdir, file)
        if os.path.isdir(d):
            dirs.append(d)
    temp=sorted(dirs, key=lambda x: os.path.getctime(x), reverse=True)[:1][0]
    return temp[-12:]

def newest(path):
    files = os.listdir(path)
    paths = [os.path.join(path, basename) for basename in files]
    return max(paths, key=os.path.getctime)

if __name__=="__main__":
    try:
        config,restart = sys.argv[1:]
        config = [str(config)]
    except ValueError:
        pass
        try:
            config = sys.argv[1:]
            restart=0
        except :
            print('Something went wrong')
    fibers_OP(config,restart)
     