import sys 
import xarray as xr
import numpy as np
import os 
import math
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

coords = xr.open_dataset('/home/jvalenti/MOAD/grid/coordinates_seagrid_SalishSea201702.nc', decode_times=False) #Here we get the grid
mask = xr.open_dataset('/home/jvalenti/MOAD/grid2/mesh_mask202108_TDV.nc') #Here we get the depth levels


def count_inside_grid_cell(center_x, center_y, cell_width, cell_height, x,y,Ni):
    deg2met = 111319.5
    cell_width = cell_width/(deg2met*np.cos(center_y*(math.pi/180)))
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

def conc_OP(filename,Ni=1,skipdays = 20, odt = 6, MFc = 5e6):
    '''Ni indicates over how many model grid cells we are calculating the concentrations
    odt is the frequency for output in the OP result file
    MFc is the number of Particles considered in the super-individual
    skipdays trims so many days out of the concentration calculation'''
    
    print(filename[0].split('.')[0]+'.npy')
    Ni =  int(Ni)
    ds = xr.open_dataset(filename[0],decode_times=False)
    zlevels = mask.gdepw_0[0,:,1,1].values
    DS=ds.to_dataframe()
    DS = DS[DS.status==1]
    DS = DS[DS.time>skipdays*86400] #skip 20 first days of run
    tmax = ((np.max(DS.time)-(skipdays*86400))/odt*3600) #output timesteps left
    print(f'Calculating average over {tmax} timesteps')
    lat = coords.nav_lat[::Ni,::Ni]
    lon = coords.nav_lon[::Ni,::Ni]
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
            conc[k,:,j]+= count_inside_grid_cell(lon[:,j], lat[:,j], cell_width[:,j], cell_height[:,j],X,Y,Ni)*MFc/(BOXarea*tmax) #Averaged over remaining timesteps
    for k in range(len(zlevels)-1):
        conc[k,:,:] = (conc[k,:,:]-conc[k+1,:,:])/(zlevels[k+1]-zlevels[k]) #Calculate concentration at each depth layer
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
     
