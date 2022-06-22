import sys
import xarray as xr
import numpy as np
import os
import yaml
import math
from datetime import datetime, timedelta
from parcels import FieldSet, Field, VectorField, ParticleSet, JITParticle, ErrorCode, ParcelsRandom, Variable

sys.path.append('/home/jvalenti/MOAD/analysis-jose/Source')
from Kernels_1 import DeleteParticle, Buoyancy, AdvectionRK4_3D, Stokes_drift, Beaching, Unbeaching, turb_mix, Biofilm
from OP_functions import *

def fibers_OP(config,local=0,restart=0):
    param = load_config(config)
    #Definitions
    start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
    length = param['param']['length'] # Set Time length [days] 
    dt = param['param']['dt'] #toggle between - or + to pick backwards or forwards 
    N = param['param']['N'] # number of deploying locations
    n = param['param']['n'] # 1000   # number of particles per location
    dmin = param['param']['dmin'] #minimum depth
    dd = param['param']['dd'] #max depth difference from dmin
    name = param['file']['name'] #name output file
    dtp = param['param']['dtp'] #how often particle released in hours
    rrr = param['param']['r'] #radious of particle deployment
# Define paths
    paths = path(local)
#Set outfall coordinates (Modify to choose other deploying location)    
    coord=xr.open_dataset(paths['coords'],decode_times=False)
    #outf_lat=coord['nav_lat'][445,304]
    #outf_lon=coord['nav_lon'][445,304]
    clat,clon = [49.195796], [-122.913127]
    #clon, clat = [float(outf_lon)],[float(outf_lat)] 
    duration = timedelta(days=length)
    x_offset, y_offset, z = p_deploy(N,n,dmin,dd,rrr)

#Set deploy locations
    lon = np.zeros([N,n])
    lat = np.zeros([N,n])
    for i in range(N):
        lon[i,:]=(clon[i] + x_offset[i,:])
        lat[i,:]=(clat[i] + y_offset[i,:])

#Set start date time and the name of the output file

    daterange = [start+timedelta(days=i) for i in range(length)]
    fn =  name + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.nc'
    outfile = os.path.join(paths['out'], fn)

####BUILD FIELDS FOR SIMULATION######

#Fill in the list of variables that you want to use as fields
    varlist=['U','V','W','R','T']
    filenames,variables,dimensions=filename_set(start,length,varlist,local)
    field_set=FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=True)

    varlist=['US','VS','WL']
    filenames,variables,dimensions=filename_set(start,length,varlist,local)

    us = Field.from_netcdf(filenames['US'], variables['US'], dimensions,allow_time_extrapolation=True)
    vs = Field.from_netcdf(filenames['VS'], variables['VS'], dimensions,allow_time_extrapolation=True)
    wl = Field.from_netcdf(filenames['WL'], variables['WL'], dimensions,allow_time_extrapolation=True)
    field_set.add_field(us)
    field_set.add_field(vs)
    field_set.add_field(wl)
    field_set.add_vector_field(VectorField("stokes", us, vs, wl))

    filenames,variables,dimensions=filename_set(start,length,['Bathy'],local)
    Bth = Field.from_netcdf(filenames['Bathy'], variables['Bathy'], dimensions,allow_time_extrapolation=True)
    field_set.add_field(Bth)

    filenames,variables,dimensions=filename_set(start,length,['FS'],local)
    Fraser = Field.from_netcdf(filenames['FS'], variables['FS'], dimensions,allow_time_extrapolation=True,timestamps=get_timestamps(start,length))
    field_set.add_field(Fraser)

    filenames,variables,dimensions=filename_set(start,length,['Kz'],local)
    Kz = Field.from_netcdf(filenames['Kz'], variables['Kz'], dimensions,allow_time_extrapolation=True)
    field_set.add_field(Kz)

    varlist=['Diat','Flag']
    filenames,variables,dimensions=filename_set(start,length,varlist,local)
    #MZ = Field.from_netcdf(filenames['MZ'], variables['MZ'], dimensions,allow_time_extrapolation=True)
    Diat = Field.from_netcdf(filenames['Diat'], variables['Diat'], dimensions,allow_time_extrapolation=True)
    Flag = Field.from_netcdf(filenames['Flag'], variables['Flag'], dimensions,allow_time_extrapolation=True)
    #field_set.add_field(MZ)
    field_set.add_field(Diat)
    field_set.add_field(Flag)

    MPParticle = particle_maker(param)

    
######RUN OCEAN PARCELS WITH DEFINED PARTICLE AND PRESET FIELDS
    if restart==1:
        name_temp=find_temp(paths['out'])
        os.system(f"cd {paths['out']} && parcels_convert_npydir_to_netcdf {name_temp}")
        outfile=newest(paths['out'])
        
        pset = ParticleSet.from_particlefile(field_set, MPParticle,outfile)
    else:
        pset = ParticleSet.from_list(field_set, MPParticle, lon=lon, lat=lat, depth=z, repeatdt = timedelta(hours=dtp))
    

    KERNELS = kernel_asem(pset,param)
    pset.execute(KERNELS,
                runtime=duration, 
                dt=dt,
                output_file=pset.ParticleFile(name=outfile, outputdt=timedelta(hours=1)),
                recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

def load_config(config_yaml):
   with open(config_yaml[0]) as f:
       config = yaml.safe_load(f)
   return config

def kernel_asem(pset,config):
    KER = AdvectionRK4_3D 
    if 'Buoyancy' in config['kernel']:
        KER += pset.Kernel(Buoyancy)
    if 'Stokes_drift' in config['kernel']:
        KER += pset.Kernel(Stokes_drift)
    if 'Beaching' in config['kernel']:
        KER += pset.Kernel(Beaching)
        KER += pset.Kernel(Unbeaching)
    if 'Turb_mix' in config['kernel']:
        KER += pset.Kernel(turb_mix)
    if 'Biofilm' in config['kernel']:
        KER += pset.Kernel(Biofilm)
    return KER

def particle_maker(config):
    #Define particle properties 
    class MPParticle(JITParticle):
        if 'ro' in config['particle']:        
            ro = Variable('ro', initial = config['particle']['ro'])  # config['particle']['ro']
        if 'diameter' in config['particle']:           
            diameter = Variable('diameter', initial = config['particle']['diameter'])
        if 'SDD' in config['particle']:  
            Sdd = Variable('SDD', initial = config['particle']['SDD'])
        if 'SDL' in config['particle']:           
            Sdl = Variable('SDL', initial = config['particle']['SDL'])
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
            tau = Variable('tau', initial =  config['particle']['tau']) # track number of particles
        if 'fratio' in config['particle']:  
            fratio = Variable('fratio', initial =  config['particle']['fratio']) # track number of particles
        if 'Nbac' in config['particle']:  
            Nbac = Variable('Nbac', initial =  0) # number of bacteria attached
        if 'Nflag' in config['particle']:  
            Nflag = Variable('Nflag', initial =  0) # number of flagellates grazing on the attached bacteria
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
     