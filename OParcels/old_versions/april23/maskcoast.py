import sys
import os
import numpy as np
import pandas as pd
import xarray as xr
sys.path.append('/home/jvalenti/MOAD/analysis-jose/Source')
from OP_functions import *

def maskcoast(local=0):
    local = 0 #Set to 0 when working on server
    paths = path(local)

    mask = xr.open_dataset(paths['mask'])

    cmask = np.zeros_like(mask.umask)
    for k in range(cmask.shape[1]):
        print(f'{100*(k/cmask.shape[1])}% completed')
        for j in range(cmask.shape[2]):
            for i in range(cmask.shape[3]):
                if mask.umask[0,k,j,i]==1:
                    if mask.umask[0,k,j-1,i-1]==0 or mask.umask[0,k,j-1,i]==0 or mask.umask[0,k,j-1,i+1]==0 or mask.umask[0,k,j,i-1]==0 or mask.umask[0,k,j,i+1]==0 or mask.umask[0,k,j+1,i-1]==0 or mask.umask[0,k,j+1,i]==0 or mask.umask[0,k,j+1,i+1]==0:
                        cmask[0,k,j,i]=1
    mask['coast_mask']=(['t','z','y','x'],cmask)
    mask.to_netcdf('/ocean/jvalenti/MOAD/grid/grid/mesh_maskBatCoast201702.nc')


if __name__=="__main__":
    try:
        config,restart = sys.argv[1:]
        config = [str(config)]
    except ValueError:
        pass
        try:
            config = sys.argv[1:]
            restart=0
        except :
            print('Something went wrong')
    maskcoast()
     