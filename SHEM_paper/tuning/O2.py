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
        if row.Depth >= 0 and j>0:
            diff = mask.gdept_0[0,:,j,i].values - row.Depth
            dd.append(diff[diff<0].argmax())
            za.append(mask.gdept_0[0,:,j,i].values[dd[-1]])
            zb.append(mask.gdept_0[0,:,j,i].values[dd[-1]+1])
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
        return N_shallow + (N_deep - N_shallow) * (z_obs - z_shallow) / (z_deep - z_shallow)

    path = '/home/jvalenti/scratch/run_SHEM/tuning/'+config[0]+'/'
    N_model = []
    for row in df.itertuples(index=False):
        try:
            fn = make_filename(path,row.folder_day)
        except FileNotFoundError:
            N_model.append('NaN')
            continue
        if row.k_above == 'NaN':
            N_model.append('NaN')
        else:
            fl = xr.load_dataset(fn)
            ab = fl['dissolved_oxygen'].isel(time_counter=0,deptht=row.k_above,y=row.j,x=row.i).item()
            be = fl['dissolved_oxygen'].isel(time_counter=0,deptht=row.k_above+1,y=row.j,x=row.i).item()
            N_model.append(interp_depth(ab, be, row.z_above, row.z_bellow, row.Depth))
            fl.close()
        print(N_model[-1])
    df['O2_model'] = N_model


    df.to_excel('DO_model_'+config[0]+'.xlsx')

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print(sys.argv)
        print('Something went wrong')
    O2()