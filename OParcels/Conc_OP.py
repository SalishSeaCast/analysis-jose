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

def get_conc(latmin,latmax,lonmin,lonmax, conc):
    mask = xr.open_dataset(paths['mask'])
    z = mask.gdepw_0[0,:,240,340]
    vol=xr.open_dataset('/home/jvalenti/MOAD/grid/grid/mesh_maskBV201702.nc')['volume_cell']
    jjii = xr.open_dataset('~/MOAD/grid/grid/grid_from_lat_lon_mask999.nc')
    j = [jjii.jj.sel(lats=latmin, lons=lonmin, method='nearest').item()]
    i = [jjii.ii.sel(lats=latmin, lons=lonmin, method='nearest').item()]
    j.append(jjii.jj.sel(lats=latmin, lons=lonmax, method='nearest').item())
    i.append(jjii.ii.sel(lats=latmin, lons=lonmax, method='nearest').item())
    j.append(jjii.jj.sel(lats=latmax, lons=lonmin, method='nearest').item())
    i.append(jjii.ii.sel(lats=latmax, lons=lonmin, method='nearest').item())
    j.append(jjii.jj.sel(lats=latmax, lons=lonmax, method='nearest').item())
    i.append(jjii.ii.sel(lats=latmax, lons=lonmax, method='nearest').item())
    a=[min(j),max(j),min(i),max(i)]
    Len = (a[1]-a[0])*(a[3]-a[2])
    SD = []
    Mean = []
    for ki in range(len(z)):
        values = []
        vols = []
        for j in range(a[0],a[1],1):
            for i in range(a[2],a[3],1):
                values.append(conc[j,i,ki])
                vols.append(vol[0,ki,j,i])
        values = np.array(values)
        vols = np.array(vols)
        Mean.append(np.sum(values)/np.sum(vols))
        valuess = np.divide(values,vols)
        SD.append(np.std(valuess)/np.sqrt(Len))
    return Mean,SD
def MP_measure(conc):
    print('0% done')
    PugC_MP,PugC_SE=get_conc(47.65,47.75,-122.5,-122.3,conc)
    print('15% done')
    PugN_MP,PugN_SE=get_conc(48,48.1,-122.75,-122.55,conc)
    print('30% done')
    JdFE_MP,JdFE_SE=get_conc(48.2,48.3,-123.25,-123.05,conc)
    print('45% done')
    JdFW_MP,JdFE_SW=get_conc(48.3,48.4,-124.25,-124.05,conc)
    print('60% done')
    SoGC_MP,SoGC_SE =get_conc(49.3,49.4,-124,-123.8,conc)
    print('75% done')
    SoGN_MP,SoGN_SE =get_conc(49.8,49.9,-124.9,-124.7,conc)
    print('90% done')
    Fraser_MP,Fraser_SE =get_conc(49,49.1,-123.5,-123.3,conc)
    print('100% done')
    return PugC_MP,PugC_SE,PugN_MP,PugN_SE,JdFE_MP,JdFE_SE,JdFW_MP,JdFE_SW,SoGC_MP,SoGC_SE,SoGN_MP,SoGN_SE,Fraser_MP,Fraser_SE


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
    MFc = param['param']['MFc']
    return outfile, MFc

def Conc_OP(config):
    outfile, MFc=loadyamls(config)
    coords=xr.open_dataset(paths['coords'],decode_times=False)
    ds = xr.open_dataset(outfile)
    DS=ds.to_dataframe()
    time=np.array(DS.xs(0, level='traj').iloc[:,3])
    mask = xr.open_dataset(paths['mask'])
    conc=np.zeros((coords.nav_lon.shape[0],coords.nav_lon.shape[1],mask.gdepw_0.shape[1]))
    jjii = xr.open_dataset('~/MOAD/grid/grid/grid_from_lat_lon_mask999.nc')
    arr = mask.gdepw_0[0,:,240,340]
    dss=DS[DS.beached==0]## In the water column
    dssla=np.array(dss.lat)
    dsslo=np.array(dss.lon)
    dsscon= float(MFc)
    dssdep=np.array(dss.z)
    
    # for i in range(len(dss)):
    #     if i%5000==0:
    #         print(f'{100*i/len(dss):.2f}% done.')
    #     jj = jjii.jj.sel(lats=dssla[i], lons=dsslo[i], method='nearest').item()
    #     ii = jjii.ii.sel(lats=dssla[i], lons=dsslo[i], method='nearest').item()
    #     try:
    #         dep = (np.abs(arr - dssdep[i])).argmin()
    #         if arr[dep] > dssdep[i]:
    #             dep+=-1
    #         conc[jj,ii,dep] += dsscon
    #     except ValueError:
    #         pass
        
    # data_set=xr.Dataset(coords={'lat': (['x', 'y'], coords.nav_lat.data),
    #                 'lon': (['x', 'y'], coords.nav_lon.data),'depth':arr})
    # data_set["Prob"]=(['x', 'y','z'], conc)
    conc=xr.open_dataarray('/home/jvalenti/MOAD/results/Salish2018_long_prob2018.nc')
    #param = load_config(config)
    #data_set.load().to_netcdf(path='/home/jvalenti/MOAD/results/'+param['file']['name']+'_prob'+str(param['startdate']['year'])+'.nc')
    PugC_MP,PugC_SE,PugN_MP,PugN_SE,JdFE_MP,JdFE_SE,JdFW_MP,JdFE_SW,SoGC_MP,SoGC_SE,SoGN_MP,SoGN_SE,Fraser_MP,Fraser_SE=MP_measure(conc)
    dict = {'PugC_MP':PugC_MP,'PugC_SE':PugC_SE,'PugN_MP':PugN_MP,'PugN_SE':PugN_SE,'JdFE_MP':JdFE_MP,'JdFE_SE':JdFE_SE,'JdFW_MP':JdFW_MP,
    'JdFE_SW':JdFE_SW,'SoGC_MP':SoGC_MP,'SoGC_SE':SoGC_SE,'SoGN_MP':SoGN_MP,'SoGN_SE':SoGN_SE,'Fraser_MP':Fraser_MP,'Fraser_SE':Fraser_SE}
    df = pd.DataFrame(dict) 
    df.to_csv('resultsSalish2.csv')
    

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print('Something went wrong')
    Conc_OP(config)