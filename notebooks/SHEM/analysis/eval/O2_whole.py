import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import xarray as xr
import os
from datetime import datetime, timedelta
import sys


def O2_whole():
    df = pd.read_csv('/home/jvalenti/MOAD/analysis-jose/notebooks/SHEM/eval/dfo-2023-2024-oxygen.csv',parse_dates=['time'],index_col=False)
    df = df[df['DOXYZZ01']>0]
    df = df[df['depth']>0]

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
    
    #Make it easy to check values in the model finding the box and folder
    df['folder_day'] = df['time'].dt.strftime('%d%b%y').str.lower()
    jj = []
    ii = []
    dd = []
    za = []
    zb = []
    for longitude, group in df.groupby('longitude'):
        for latitude,subgroup in group.groupby('latitude'):
            j,i = finder(latitude,longitude)
            if j>0: #Is it in the Domain?\
                for idx, row in subgroup.iterrows():
                    jj.append(j)
                    ii.append(i)
                    if row.depth >= 0.51: #Is it below the first T level?
                        diff = mask.gdept_0[0,:,j,i].values - row.depth
                        dd.append(diff[diff<=0].argmax())
                        if int(dd[-1]) >= 39: #Is it below the last T level?
                            print('depth is deeper than model bottom')
                            za.append(mask.gdept_0[0,:,j,i].values[dd[-1]])
                            zb.append(np.nan)
                        else:
                            za.append(mask.gdept_0[0,:,j,i].values[dd[-1]])
                            zb.append(mask.gdept_0[0,:,j,i].values[dd[-1]+1])
                            tz = mask.tmask[0,dd[-1]:dd[-1]+2,j,i]
                            if tz[0] == 0:
                                za[-1] = np.nan
                            if tz[1] == 0:
                                zb[-1] = np.nan
                    else: #Choose surface T level
                        dd.append(0.5)
                        za.append(0.5)
                        zb.append(1) 
            else: #Any value works here, we are deleting these.
                for idx, row in subgroup.iterrows():
                    jj.append(j)
                    ii.append(i)
                    dd.append(0)
                    za.append(0)
                    zb.append(0) 
    

    df['j'] = jj
    df['i'] = ii
    df['k_above'] = dd
    df['z_above'] = za
    df['z_bellow'] = zb

    df = df[~np.isnan(df['z_above'])].reset_index(drop=True)
    df = df[df['j']>0].reset_index(drop=True)
    df['k_above'] = df['k_above'].astype(int)

    def interp_depth(N_shallow, N_deep, z_shallow, z_deep, z_obs):
        if N_deep > 0:
            return N_shallow + (N_deep - N_shallow) * (z_obs - z_shallow) / (z_deep - z_shallow)
        else:
            return N_shallow
    
    path = f'/home/jvalenti/scratch/run_SHEM/long_run'
    N_model = np.full(len(df), np.nan)
    print(df)
    for folder_day, group in df.groupby('folder_day'):
        print(folder_day)
        try:
            fn = make_filename(path,folder_day,res = 'h')
        except FileNotFoundError:
            continue
        with xr.open_dataset(fn, engine='h5netcdf') as ds:
            var = ds['dissolved_oxygen'].isel(time_counter=0)
            for idx, row in group.iterrows():
                ab = var.isel(deptht=slice(row.k_above, row.k_above + 2),y=row.j,x=row.i).values
                try:
                    N_model[idx] = interp_depth(ab[0], ab[1],row.z_above, row.z_bellow,row.depth)
                except IndexError:
                    print('Deeper than domain get last level')
                    N_model[idx] = ab[0]

        print(folder_day)
    df['DO_model'] = N_model
    df.to_excel('DO_model_dfo.xlsx')

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print(sys.argv)
        print('Something went wrong')
    O2_whole()