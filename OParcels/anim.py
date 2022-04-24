import sys
import os
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt, animation, rc
from datetime import datetime, timedelta
import yaml
sys.path.append('/home/jvalenti/MOAD/analysis-jose/notebooks/parcels')
from OP_functions_biofilm import *

def anim(config,fps,local=0,):
    param = load_config(config)
    #Definitions
    start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day']) #Start date
    Tmax = param['param']['length'] # Set Time length [days] 
    duration = timedelta(days=Tmax)
    dt = param['param']['dt'] #toggle between - or + to pick backwards or forwards 
    N = param['param']['N'] # number of deploying locations
    n = param['param']['n'] # 1000   # number of particles per location
    dmin = param['param']['dmin'] #minimum depth
    dd = param['param']['dd'] #max depth difference from dmin
    name = param['file']['name'] #name output file
    # Define paths
    paths = path(local)

    coord=xr.open_dataset(paths['coords'],decode_times=False)
    outf_lat=coord['nav_lat'][445,304]
    outf_lon=coord['nav_lon'][445,304]
    clon, clat = [float(outf_lon)],[float(outf_lat)] 

    daterange = [start+timedelta(days=i) for i in range(Tmax)]
    fn =  name + '_'.join(d.strftime('%Y%m%d')+'_1n' for d in [start, start+duration]) + '.nc'
    outfile = os.path.join(paths['out'], fn)


    coords = xr.open_dataset(paths['coords'], decode_times=False)
    mask = xr.open_dataset(paths['mask'])
    ds = xr.open_dataset(outfile)

    anim = mapanimationd(outfile,N,n,clon,clat,fps,local)
    f = r"/home/jvalenti/MOAD/animations/"+name+".gif" 
    FFwriter = animation.FFMpegWriter()
    anim.save(f, writer = FFwriter)

def load_config(config_yaml):
   with open(config_yaml[0]) as f:
       config = yaml.safe_load(f)
   return config

if __name__=="__main__":
    try:
        config,fps = sys.argv[1:]
        config = [str(config)]
    except ValueError:
        pass
        try:
            config = sys.argv[1:]
            fps=240
        except :
            print('Something went wrong')
    anim(config,int(fps))