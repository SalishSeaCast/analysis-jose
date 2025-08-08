import sys 
import xarray as xr
import numpy as np
import os 
import math

coords = xr.open_dataset('/home/jvalenti/MOAD/grid/coordinates_seagrid_SalishSea201702.nc', decode_times=False)
mask = xr.open_dataset('/home/jvalenti/MOAD/grid2/mesh_mask202108_TDV.nc')

def latT(lat):
    return np.cos(lat*(math.pi/180))

def count_inside_grid_cell(center_x, center_y, cell_width, cell_height, x,y,Ni):
    deg2met = 111319.5
    cell_width = cell_width/(deg2met*latT(center_y))
    cell_height = cell_height/deg2met
    min_x = np.array(center_x - cell_width / 2)
    min_y = np.array(center_y - cell_height / 2)
    max_y = np.array(center_y + cell_height*(Ni-0.5))
    max_x = np.array(center_x + cell_width*(Ni-0.5))
    
    
    inside_mask = np.logical_and.reduce([
        min_x[:, np.newaxis] <= x,
        x <= max_x[:, np.newaxis],
        min_y[:, np.newaxis] <= y,
        y <= max_y[:, np.newaxis]
    ])
    c = np.sum(inside_mask,axis=1)
    return c

def conc_OP(filename,Ni=1):
    print(filename[0].split('.')[0]+'.npy')
    Ni =  int(Ni)
    ds = xr.open_dataset(filename[0],decode_times=False)
    zlevels = mask.gdepw_0[0,:,1,1].values
    DS=ds.to_dataframe()
    DS = DS[DS.status==1]
    DS = DS[DS.time>1728000] #last 10 days of run
    tmax = int((np.max(DS.time)-1728000)/21600)
    print(f'Calculating average over {tmax} timesteps')
    MFc = 5e6/tmax #Averaged over remaining days
    lat = coords.nav_lat[::Ni,::Ni]
    lon = coords.nav_lon[::Ni,::Ni]
    td = mask.totaldepth
    cell_width = coords.e1t[0,::Ni,::Ni]
    cell_height = coords.e2t[0,::Ni,::Ni]
    x = np.array(DS.lon)
    y = np.array(DS.lat)
    z = np.array(DS.z)

    conc = np.zeros([len(zlevels),len(coords.nav_lon[::Ni,0]),len(coords.nav_lon[0,::Ni])])
    for k in range(len(zlevels)):
        print(f'{k} level starting.')  
        for j in range(len(coords.nav_lon[0,::Ni])):
            zmin = int(zlevels[k])
            X = x[z >= zmin]
            Y = y[z >= zmin]
            BOXarea = (cell_width[:,j]* cell_height[:,j])*Ni**2
            #BOXvolume = (cell_width[j,:]* cell_height[j,:]*(td[j,:]-zmin))
            conc[k,:,j]+= count_inside_grid_cell(lon[:,j], lat[:,j], cell_width[:,j], cell_height[:,j],X,Y,Ni)*MFc/BOXarea
    for k in range(len(zlevels)-1):
        conc[k,:,:] = (conc[k,:,:]-conc[k+1,:,:])/(zlevels[k+1]-zlevels[k])
    np.save(filename[0].split('.')[0]+'.npy',conc)

if __name__=="__main__":
    try:
        filename,Ni = sys.argv[1:]
        filename = [str(filename)]
    except ValueError:
        print('Not reducing')
        try:
            filename = sys.argv[1:]
            restart=0
        except :
            print('Something went wrong')
    conc_OP(filename,Ni)
     
