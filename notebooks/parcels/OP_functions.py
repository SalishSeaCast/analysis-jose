import numpy as np
import os
from cartopy import crs, feature
from matplotlib import pyplot as plt, animation, rc
import xarray as xr
from datetime import datetime, timedelta

def make_prefix(date, path, res='h'):
    """Construct path prefix for local SalishSeaCast results given date object and paths dict
    e.g., /results2/SalishSea/nowcast-green.201905/daymonthyear/SalishSea_1h_yyyymmdd_yyyymmdd
    """

    datestr = '_'.join(np.repeat(date.strftime('%Y%m%d'), 2))
    folder = date.strftime("%d%b%y").lower()
    prefix = os.path.join(path, f'{folder}/SalishSea_1{res}_{datestr}')
    
    return prefix
colores=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def scatter_particles(ax, N ,n,nmin, nmax,yvar,lon,HD=0,colors=colores):
    scatter=[]
    #N is the number of stations, n the number of particles
    #Color is a list of strings picking the desired colors
    #nmin is t0, nmax is tmax
    #yvar is the y coordinate, lon is the longitud array
    #HD (1 for cartopy plot) 0 otherwise
    starts = np.arange(0,N*n,n)
    ends = np.arange(n-1,N*n,n)
    if N < len(colors):
        colors = colors[0:N]
    elif N > len(colors):
        con = 0
        while N > len(colors):
            colors.append(colors[con])
            con+=1
    if HD == 0:
        if nmin==nmax:
            for i in range(N):
                scatter.append(ax.scatter(lon[starts[i]:ends[i], nmin], yvar[starts[i]:ends[i], nmin],c=colors[i],s=5))
        else:
            for i in range(N):
                scatter.append(ax.scatter(lon[starts[i]:ends[i], nmin:nmax], yvar[starts[i]:ends[i], nmin:nmax],c=colors[i],s=5))
    else:
        if nmin==nmax:
            for i in range(N):
                scatter.append(ax.scatter(lon[starts[i]:ends[i], nmin], yvar[starts[i]:ends[i], nmin],c=colors[i],s=5,transform=crs.PlateCarree(),zorder=2))
        else:
            for i in range(N):
                scatter.append(ax.scatter(lon[starts[i]:ends[i], nmin:nmax], yvar[starts[i]:ends[i], nmin:nmax],c=colors[i],s=5,transform=crs.PlateCarree(),zorder=2))
        
    return scatter

ss=[]

def mapanimation(outfile,N,n,clon,clat,dmin,dmax,fps=1):
    coords,mask,ds=output(outfile)
    fig = plt.figure(figsize=(19, 8))
    ax = plt.axes(xlim=(-127,-121),ylim=(46.8,51.2))
    ax.contour(coords.nav_lon, coords.nav_lat, mask.mbathy[0,:,:],colors='k',linewidths=0.1)
    ax.contourf(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='lightgray')
    ax.contour(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='k')
    ax.grid()
    ax.set_aspect(1/1)
    plt.ylabel('Latitude',fontsize=16)
    plt.xlabel('Longitude',fontsize=16)
    t = ax.text(0.02, 0.02, '', transform=ax.transAxes)
    t.set_text('')
    ss = scatter_particles(ax, N,n, 0,0, ds.lat,ds.lon)

    def update(frame):
        tstamp = ds.time[0, frame].values.astype('datetime64[s]').astype(datetime)
        t.set_text(tstamp.strftime('%Y-%b-%d %H:%M UTC'))
        global ss
        for scat in ss:
            scat.remove()
        ss = scatter_particles(ax, N,n, frame,frame, ds.lat,ds.lon)
        ss.append(ax.scatter(clon,clat,c='r', marker='*', linewidths=2))
        return ss

    anim= animation.FuncAnimation(fig, update, frames=np.arange(0,len(ds.lon[0,:]),fps))
    return anim

def visual(outfile,N,n,clon,clat,dmin,dmax, nmin=0, nmax=-1):
    coords,mask,ds = output(outfile)
    paths = path()
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(19, 8))
    ax1.contour(coords.nav_lon, coords.nav_lat, mask.mbathy[0,:,:],colors='k',linewidths=0.1)
    ax1.contourf(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='lightgray')
    ax1.contour(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='k')

    scatter_particles(ax1, N,n, nmin, nmax, ds.lat,ds.lon)
    ax1.scatter(clon,clat,c='g', marker='*', linewidths=1)

    scatter_particles(ax2, N,n, nmin, nmax, -ds.z,ds.lon)
    ax2.grid()
    plt.ylabel('Depth [m]')
    plt.xlabel('Longitude')
    zsc=np.repeat((dmin-dmax)/2,len(clon))
    ax2.scatter(clon,zsc,c='r', marker='*', linewidths=1)
    
def path(local = 1):
    if local == 1:
        path = {'NEMO': '/results2/SalishSea/nowcast-green.201905/',
        'coords': '/Users/jvalenti/MOAD/SSC_masks/coordinates_seagrid_SalishSea201702.nc',
        'mask': '/Users/jvalenti/MOAD/SSC_masks/mesh_mask201702.nc',
        'out': '/Users/jvalenti/MOAD/analysis-jose/notebooks/results/'}    
    else:
        path = {'NEMO': '/results2/SalishSea/nowcast-green.201905/',
        'coords': '/ocean/jvalenti/MOAD/grid/coordinates_seagrid_SalishSea201702.nc',
        'mask': '/ocean/jvalenti/MOAD/grid/mesh_mask201702.nc',
        'out': '/home/jvalenti/MOAD/analysis-jose/notebooks/results',
        'anim': '/home/jvalenti/MOAD/animations'}
    return path
    

def output(outfile):
    paths = path()
    coords = xr.open_dataset(paths['coords'], decode_times=False)
    mask = xr.open_dataset(paths['mask'])
    ds = xr.open_dataset(outfile)
    return  coords,mask,ds

def profile(N,n,length,outfile,levels=20):
    coords,mask,ds = output(outfile)
    Z = np.linspace(0,430,levels)
    starts = np.arange(0,N*n,n)
    ends = np.arange(n-1,N*n,n)
    time = length*24+1
    zn = np.zeros([len(Z)-1,time])
    for j in range(time):
        zn[:,j],z_levels = np.histogram(ds.z[:, j], bins=Z)
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(xlim=(-5,np.max(zn[:,0]+5)),ylim=(-500,0))
    plt.plot(zn[:,0],-z_levels[1:],'--',label='$t_0$')
    ax.grid()
    plt.ylabel('Depth [m]',fontsize=16)
    plt.xlabel('Particles',fontsize=16)
    plt.plot(zn[:,-1],-z_levels[1:],'--',label='$t_{end}$')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    
    
def filename_set(start,duration,varlist=['U','V','W']):
    #Build filenames
    paths = path()
    Tlist,Ulist, Vlist, Wlist = [], [], [], []
   
    for day in range(duration.days):
        path_NEMO = make_prefix(start + timedelta(days=day), paths['NEMO'])
        Ulist.append(path_NEMO + '_grid_U.nc')
        Vlist.append(path_NEMO + '_grid_V.nc')
        Wlist.append(path_NEMO + '_grid_W.nc')
        Tlist.append(path_NEMO + '_grid_T.nc')

    # Load NEMO forcing 
    filenames = {
        'U': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Ulist},
        'V': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Vlist},
        'W': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Wlist},
        'T': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Tlist},
        'S': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Tlist},
    }
    variables = {'U': 'vozocrtx', 'V': 'vomecrty','W': 'vovecrtz','T':'votemper','S':'vosaline'}
    dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw','time': 'time_counter'}
    file2,var2 = {},{}
    for var in varlist:
        file2[var]=filenames[var]
        var2[var]=variables[var]
    return file2,var2,dimensions

