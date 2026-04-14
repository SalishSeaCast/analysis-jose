import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import xarray as xr
import os
from datetime import datetime, timedelta
import sys


def O2():
    df = pd.read_excel('nutrients_2018dfo.xlsx',parse_dates=['dtUTC'],index_col=0)
    df = df[(df['dtUTC']>=dt.datetime(2018,2,27)) & (df['dtUTC']<=dt.datetime(2018,7,1))].reset_index(drop=True)
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
    df['folder_day'] = df['dtUTC'].dt.strftime('%d%b%y').str.lower()
    jj = []
    ii = []
    dd = []
    za = []
    zb = []
    for row in df.itertuples(index=False):
        j,i = finder(row.Lat,row.Lon)
        jj.append(j)
        ii.append(i)
        if row.Depth >= 0 and j>0 and row.Depth != 0.5:
            diff = mask.gdept_0[0,:,j,i].values - row.Depth
            dd.append(diff[diff<0].argmax())
            za.append(mask.gdept_0[0,:,j,i].values[dd[-1]])
            zb.append(mask.gdept_0[0,:,j,i].values[dd[-1]+1])
        elif row.Depth == 0.5 and j>0:
            dd.append(0)
            za.append(0)
            zb.append(1)
        else:
            dd.append('NaN')
            za.append('NaN')
            zb.append('NaN')
    df['j'] = jj
    df['i'] = ii
    df['k_above'] = dd
    df['z_above'] = za
    df['z_bellow'] = zb

    def interp_depth(N_shallow, N_deep, z_shallow, z_deep, z_obs):
        if N_deep>0:
            return N_shallow + (N_deep - N_shallow) * (z_obs - z_shallow) / (z_deep - z_shallow)
        else:
            return N_shallow

    runs = ['SSBase','SHEM18','tuning/diat_pref','tuning/exc_hbac','tuning/exc_hbac_2','tuning/growth_flag','tuning/growth_flag_2','tuning/mort_hbac','tuning/pred_flag','tuning/remin','tuning/remin2','tuning/predmine','tuning/mort_hbac_2']
    names = ['SSBase','SHEM18','diat_pref','exc_hbac','exc_hbac_2','growth_flag','growth_flag_2','mort_hbac','pred_flag','remin','remin2','predmine','mort_hbac_2']
    for i,name in enumerate(names):
        print(f'Starting: {runs[i]}')
        path = f'/home/jvalenti/scratch/run_SHEM/{runs[i]}/'
        N_model = np.full(len(df), np.nan)

        for folder_day, group in df.groupby('folder_day'):
            try:
                if name=='SSBase':
                    fn = make_filename(path,row.folder_day,'chem_T')
                else:
                    fn = make_filename(path,row.folder_day)
            except FileNotFoundError:
                continue

            with xr.open_dataset(fn, engine='h5netcdf') as ds:
                var = ds['dissolved_oxygen'].isel(time_counter=0)

                for idx, row in group.iterrows():
                    if row.k_above == 'NaN':
                        continue
                    ab = var.isel(deptht=slice(row.k_above, row.k_above + 2),y=row.j,x=row.i).values
                    N_model[idx] = interp_depth(ab[0], ab[1],row.z_above, row.z_bellow,row.Depth)
            print(folder_day)
        df[name] = N_model
    df.to_excel('DO_model_whole.xlsx')

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print(sys.argv)
        print('Something went wrong')
    O2()