import sys
import os
import numpy as np
import pandas as pd
import xarray as xr

import yaml

sys.path.append('/home/jvalenti/MOAD/analysis-jose/Source')
from OP_functions import *

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
    daterange = [start+timedelta(days=i) for i in range(Tmax)]
    fn =  name + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.nc'
    outfile = os.path.join(paths['out'], fn)
    return outfile

def Prob_OP(config):
    outfile=loadyamls(config)
    coords=xr.open_dataset(paths['coords'],decode_times=False)
    ds = xr.open_dataset(outfile)
    DS=ds.to_dataframe()
    time=np.array(DS.xs(0, level='traj').iloc[:,1])
    mask = xr.open_dataset(paths['mask'])
    conc=np.zeros((coords.nav_lon.shape[0],coords.nav_lon.shape[1],mask.gdepw_0.shape[1]))
    jjii = xr.open_dataset('~/MOAD/grid/grid/grid_from_lat_lon_mask999.nc')
    arr = mask.gdepw_0[0,:,240,340]
    dss=DS[DS.beached==0]## In the water column
    dssla=np.array(dss.lat)
    dsslo=np.array(dss.lon)
    dsscon= np.array(dss.tau)
    dssdep=np.array(dss.z)
    
    for i in range(len(dss)):
        if i%5000==0:
            print(f'{100*i/len(dss):.2f}% done.')
        jj = jjii.jj.sel(lats=dssla[i], lons=dsslo[i], method='nearest').item()
        ii = jjii.ii.sel(lats=dssla[i], lons=dsslo[i], method='nearest').item()
        try:
            dep = (np.abs(arr - dssdep[i])).argmin()
            if arr[dep] > dssdep[i]:
                dep+=-1
            conc[jj,ii,dep] += dsscon[i]
        except ValueError:
            pass
        
    data_set=xr.Dataset(coords={'lat': (['x', 'y'], coords.nav_lat.data),
                    'lon': (['x', 'y'], coords.nav_lon.data),'depth':arr})
    data_set["Prob"]=(['x', 'y','z'], conc)
    param = load_config(config)
    data_set.load().to_netcdf(path='/home/jvalenti/MOAD/results/'+param['file']['name']+'_prob'+str(param['startdate']['year'])+'.nc')

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print('Something went wrong')
    Prob_OP(config)