import sys
import os
import numpy as np
import pandas as pd
import xarray as xr

import matplotlib
from matplotlib import pyplot as plt, animation, rc,colors
from datetime import datetime, timedelta
from cartopy import crs, feature
import cmocean
import yaml

from IPython.display import Image
rc('animation', html='html5')

sys.path.append('/home/jvalenti/MOAD/analysis-jose/Source')
from OP_functions import *

# Define paths
local = 0 #Set to 0 when working on server
paths = path(local)
coords = xr.open_dataset(paths['coords'], decode_times=False)
mask = xr.open_dataset(paths['mask'])

config='/home/jvalenti/MOAD/analysis-jose/OParcels/Beach.yaml'
param = load_config1(config)
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

config='/home/jvalenti/MOAD/analysis-jose/OParcels/Beach2.yaml'
param = load_config1(config)
start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
Tmax = param['param']['length'] # Set Time length [days] 
duration = timedelta(days=Tmax)
name = param['file']['name'] #name output file
fn =  name + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.nc'
outfile2 = os.path.join(paths['out'], fn)

config='/home/jvalenti/MOAD/analysis-jose/OParcels/Beach3.yaml'
param = load_config1(config)
start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
Tmax = param['param']['length'] # Set Time length [days] 
duration = timedelta(days=Tmax)
name = param['file']['name'] #name output file
fn =  name + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.nc'
outfile3 = os.path.join(paths['out'], fn)

config='/home/jvalenti/MOAD/analysis-jose/OParcels/Beach4.yaml'
param = load_config1(config)
start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
Tmax = param['param']['length'] # Set Time length [days] 
duration = timedelta(days=Tmax)
name = param['file']['name'] #name output file
fn =  name + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.nc'
outfile4 = os.path.join(paths['out'], fn)

config='/home/jvalenti/MOAD/analysis-jose/OParcels/Beach5.yaml'
param = load_config1(config)
start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
Tmax = param['param']['length'] # Set Time length [days] 
duration = timedelta(days=Tmax)
name = param['file']['name'] #name output file
fn =  name + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.nc'
outfile5 = os.path.join(paths['out'], fn)


ds = xr.open_dataset(outfile)
ds2 = xr.open_dataset(outfile2)
ds3 = xr.open_dataset(outfile3)
ds4 = xr.open_dataset(outfile4)
ds5 = xr.open_dataset(outfile5)

a=ds.to_dataframe()
a2=ds2.to_dataframe()
a3=ds3.to_dataframe()
a4=ds4.to_dataframe()
a5=ds5.to_dataframe()

obsss = np.zeros(1921)
obsss2 = np.zeros_like(obsss)
obsss3 = np.zeros_like(obsss)
obsss4 = np.zeros_like(obsss)
obsss5 = np.zeros_like(obsss)


for i in range(12000):
    df=a.loc[i,:]
    df2=a2.loc[i,:]
    df3=a3.loc[i,:]
    df4=a4.loc[i,:]
    df5=a5.loc[i,:] 
    obsss[df[df.beached!=0].first_valid_index()]+=1
    obsss2[df2[df2.beached!=0].first_valid_index()]+=1
    obsss3[df3[df3.beached!=0].first_valid_index()]+=1
    obsss4[df4[df4.beached!=0].first_valid_index()]+=1
    obsss5[df5[df5.beached!=0].first_valid_index()]+=1
    if i%100 == 0:
        print(f'{int(100*i/12000)}% completed')

d = {"b1": obsss, "b2": obsss2,"b3":obsss3,"b4":obsss4,"b5":obsss5}
df = pd.DataFrame(d)
df.to_csv('/home/jvalenti/MOAD/analysis-jose/OParcels/winter.csv')