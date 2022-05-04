import sys
import os
import numpy as np
import pandas as pd
import xarray as xr

import yaml

sys.path.append('/home/jvalenti/MOAD/analysis-jose/notebooks/parcels')
from OP_functions_biofilm import *

# Define paths
local = 0 #Set to 0 when working on server
paths = path(local)

def load_config(config_yaml):
   with open(config_yaml[0]) as f:
       config = yaml.safe_load(f)
   return config

def loadyamls(config):
    param = load_config(config)
    start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
    Tmax = param['param']['length'] # Set Time length [days] 
    duration = timedelta(days=Tmax)
    dt = param['param']['dt'] #toggle between - or + to pick backwards or forwards 
    N = param['param']['N'] # number of deploying locations
    n = param['param']['n'] # 1000   # number of particles per location
    dmin = param['param']['dmin'] #minimum depth
    dd = param['param']['dd'] #max depth difference from dmin
    name = param['file']['name'] #name output file
    WS = param['particle']['Ws']
    daterange = [start+timedelta(days=i) for i in range(Tmax)]
    fn =  name + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.nc'
    outfile = os.path.join(paths['out'], fn)
    return outfile, WS

def concen_OP(config):
    outfile,WS=loadyamls(config)
    coords=xr.open_dataset(paths['coords'],decode_times=False)
    ds = xr.open_dataset(outfile)
    conc=np.zeros((ds.time.shape[1],coords.nav_lon.shape[0],coords.nav_lon.shape[1]))
    DS=ds.to_dataframe()
    time=np.array(DS.xs(0, level='traj').iloc[:,1])
    jjii = xr.open_dataset('~/MOAD/grid/grid/grid_from_lat_lon_mask999.nc')
    for tt,t in enumerate(time):
        dss=DS[DS.time==t]
        dssla=np.array(dss.lat)
        dsslo=np.array(dss.lon)
        for i in range(len(dss)):
            jj = jjii.jj.sel(lats=dssla[i], lons=dsslo[i], method='nearest').item()
            ii = jjii.ii.sel(lats=dssla[i], lons=dsslo[i], method='nearest').item()
            conc[tt,jj,ii] += 1
    time=np.datetime_as_string(time, unit='m')
    data_set=xr.Dataset( coords={'time': time, 'lat': (['x', 'y'], coords.nav_lat.data),
                    'lon': (['x', 'y'], coords.nav_lon.data)})
    data_set["conc"]=([ 'time','x', 'y'],  conc)
    data_set.to_netcdf(path=config[0][:-5]+'.nc')

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print('Something went wrong')
    concen_OP(config)