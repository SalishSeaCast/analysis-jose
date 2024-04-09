from glob import glob
from os import path
import xarray as xr
import sys
import zarr

def zarr_tonet(config):
    config =  str(config[0])
    name =  config.split('1n')[0]
    fname = '/scratch/jvalenti/OParcels_runs/Parcels_alpha/results/'+config
    print(fname)
    files = glob(path.join(fname, "proc*"))
    ds = xr.concat(
        [xr.open_zarr(f) for f in files],
        dim="trajectory",
        compat="no_conflicts",
        coords="minimal",
    )
    print('done')
   #ds = ds.drop_vars(['alpha','Ub','cellvol','diameter','fratio','length','tau','wa','ws','wm'])
    ds.to_netcdf('/home/jvalenti/projects/rrg-allen/jvalenti/'+name+'.nc')

if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print('Something did go wrong')
    zarr_tonet(config)
