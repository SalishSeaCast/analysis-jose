import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import xarray as xr
import os
from datetime import datetime, timedelta
import sys


def O2():
    df = pd.read_csv('PugetSound_2018.csv',parse_dates=['Sample_Date'])
    df = df[(df['Sample_Date']>=dt.datetime(2018,2,27)) & (df['Sample_Date']<=dt.datetime(2018,7,1))].reset_index(drop=True)
    df = df[df['DO']>0]

    jjii = xr.open_dataset('~/MOAD/grid/grid_from_lat_lon_mask999.nc')
    def finder(lati,loni):
        j = [jjii.jj.sel(lats=lati, lons=loni, method='nearest').item()][0]
        i = [jjii.ii.sel(lats=lati, lons=loni, method='nearest').item()][0]
        return j,i

    def make_filename(path_run,folder, var='biol_T', res='d'):
        """Construct path prefix for local SHEM results given date object and paths dict"""
        prefix = os.path.join(path_run, f'{folder}/')
        fname = []
        for file in os.listdir(prefix):
            if (var in file) and ('_1'+res) in file:
                fname.append(file)
        if len(fname)>1:
            print('more than one file found') 
        
        return os.path.join(prefix, fname[0])   

    mask = xr.open_dataset('/home/jvalenti/MOAD/grid2/mesh_mask202108_TDV.nc') 


    df['folder_day'] = df['Sample_Date'].dt.strftime('%d%b%y').str.lower()
    jj = []
    ii = []
    dd = []
    za = []
    zb = []
    for row in df.itertuples(index=False):
        j,i = finder(row.Latitude,row.Longitude)
        jj.append(j)
        ii.append(i)
        if row.Depth >= 0 and j>0:
            diff = mask.gdept_0[0,:,j,i].values - row.Depth
            dd.append(diff[diff<=1e-4].argmax())
            za.append(mask.gdept_0[0,:,j,i].values[dd[-1]])
            zb.append(mask.gdept_0[0,:,j,i].values[dd[-1]+1])
            tz = mask.tmask[0,dd[-1]:dd[-1]+2,j,i]
            if tz[0] == 0:
                za[-1] = np.nan
            if tz[1] == 0:
                zb[-1] = np.nan  
        else:
            dd.append(np.nan)
            za.append(np.nan)
            zb.append(np.nan) 
    df['j'] = jj
    df['i'] = ii
    df['k_above'] = dd
    df['z_above'] = za
    df['z_bellow'] = zb

    df = df[~np.isnan(df['z_above'])].reset_index(drop=True)
    df['k_above'] = df['k_above'].astype(int)

    def interp_depth(N_shallow, N_deep, z_shallow, z_deep, z_obs):
        if N_deep>0:
            return N_shallow + (N_deep - N_shallow) * (z_obs - z_shallow) / (z_deep - z_shallow)
        else:
            return N_shallow

    path = '/home/jvalenti/scratch/run_SHEM/tuning/'+config[0]+'/'
    if config[0] == 'SSBase' or config[0] == 'SHEM18':
        path = '/home/jvalenti/scratch/run_SHEM/'+config[0]+'/'
    N_model = np.full(len(df), np.nan)

    for folder_day, group in df.groupby('folder_day'):
        try:
            if config[0]=='SSBase':
                fn = make_filename(path,folder_day,'chem_T')
            else:
                fn = make_filename(path,folder_day)
        except FileNotFoundError:
            continue
        with xr.open_dataset(fn, engine='h5netcdf') as ds:
            var = ds['dissolved_oxygen'].isel(time_counter=0)
            for idx, row in group.iterrows():
                ab = var.isel(deptht=slice(row.k_above, row.k_above + 2),y=row.j,x=row.i).values
                N_model[idx] = interp_depth(ab[0], ab[1],row.z_above, row.z_bellow,row.Depth)
        print(folder_day)
    df['DO_model'] = N_model

    df.to_csv('DO_puget_'+config[0]+'.csv', index=False)

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print(sys.argv)
        print('Something went wrong')
    O2()