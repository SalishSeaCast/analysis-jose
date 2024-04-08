import sys 
import xarray as xr
import numpy as np
import os 
import math

coords = xr.open_dataset('/home/jvalenti/MOAD/grid/coordinates_seagrid_SalishSea201702.nc', decode_times=False)
mask = xr.open_dataset('/home/jvalenti/MOAD/grid2/mesh_mask202108_TDV.nc')

def latT(lat):
    return np.cos(lat*(math.pi/180))

def count_inside_grid_cell(center_x, center_y, cell_width, cell_height, x,y):
    deg2met = 111319.5
    cell_width = cell_width/(deg2met*latT(center_y))
    cell_height = cell_height/deg2met
    min_x = np.array(center_x - cell_width / 2)
    max_x = np.array(center_x + cell_width / 2)
    min_y = np.array(center_y - cell_height / 2)
    max_y = np.array(center_y + cell_height / 2)
    inside_mask = np.logical_and.reduce([
        min_x[:, np.newaxis] <= x,
        x <= max_x[:, np.newaxis],
        min_y[:, np.newaxis] <= y,
        y <= max_y[:, np.newaxis]
    ])
    c = np.sum(inside_mask,axis=1)
    return c



def conc_OP(filename,Zi):
    print(filename[0])
    ds = xr.open_dataset(filename[0],decode_times=False)
    MFc = 5e6
    zlevels = [0,5,10,50,100,800]
    DS=ds.to_dataframe()
    #DS = DS[DS.z < 800]
    DS = DS[DS.status==1]
    lat = coords.nav_lat
    lon = coords.nav_lon
    td = mask.totaldepth
    cell_width = coords.e1t[0,:,:]
    cell_height = coords.e2t[0,:,:]
    x = np.array(DS.lon)
    y = np.array(DS.lat)
    z = np.array(DS.z)

    conc = np.zeros((len(zlevels),coords.nav_lon.shape[0],coords.nav_lon.shape[1]))
    for k in range(len(zlevels)-2,-1,-1):
        print(f'{k} level starting.')  
        for j in range(coords.nav_lon.shape[0]):
            zmin = int(zlevels[k-1])
            zmax = int(zlevels[k])
            X = x[np.logical_and(z >= zmin, z < zmax)]
            Y = y[np.logical_and(z >= zmin, z < zmax)]
            BOXvolume = (cell_width[j,:]* cell_height[j,:]*(td[j,:]-zmin))
            conc[k,j,:]+= count_inside_grid_cell(lon[j,:], lat[j,:], cell_width[j,:], cell_height[j,:],X,Y)*MFc/BOXvolume
    np.save('concentration_31days.npy',conc)

if __name__=="__main__":
    try:
        filename,restart = sys.argv[1:]
        filename = [str(filename)]
    except ValueError:
        print('Not restarting')
        try:
            filename = sys.argv[1:]
            restart=0
        except :
            print('Something went wrong')
    conc_OP(filename)
     
