import numpy as np
import os
import sys
import math
from cartopy import crs, feature
from matplotlib import pyplot as plt, animation, rc
import xarray as xr 
import yaml
import cmocean
from datetime import datetime, timedelta
from parcels import FieldSet, Field, VectorField, ParticleSet, JITParticle, ErrorCode, ParcelsRandom, Variable

sys.path.append('/home/jvalenti/MOAD/analysis-jose/Source')
from OP_Kernels import DeleteParticle, Buoyancy, AdvectionRK4_3D, Stokes_drift, Beaching, Unbeaching, turb_mix,turb_mix2 , Biofilm
 
def path(local = 1):
    '''Change with your paths'''
    if local == 1:
        path = {'NEMO': '/Users/jvalenti/MOAD/data/',
        'coords': '/Users/jvalenti/MOAD/SSC_masks/coordinates_seagrid_SalishSea201702.nc',
        'mask': '/Users/jvalenti/MOAD/SSC_masks/mesh_mask201703.nc',
        'out': '/Users/jvalenti/MOAD/results/',
        'home': '/Users/jvalenti/MOAD/analysis-jose/notebooks/parcels',
        'anim': '/Users/jvalenti/MOAD/animations'}
    else:
        path = {'NEMO': '/results2/SalishSea/nowcast-green.201905/',
        'coords': '/ocean/jvalenti/MOAD/grid/coordinates_seagrid_SalishSea201702.nc',
        'coordsWW3': '/ocean/jvalenti/MOAD/grid/WW3_grid.nc',
        'mask': '/ocean/jvalenti/MOAD/grid/grid/mesh_maskBatCoast201702.nc',
        'out': '/home/jvalenti/MOAD/results',
        'home': '/home/jvalenti/MOAD/analysis-jose/notebooks/parcels',
        'anim': '/home/jvalenti/MOAD/animations'}
    return path


def make_prefix(date, path, res='h'):
    """Construct path prefix for local SalishSeaCast results given date object and paths dict
    e.g., /results2/SalishSea/nowcast-green.201905/daymonthyear/SalishSea_1h_yyyymmdd_yyyymmdd
    """

    datestr = '_'.join(np.repeat(date.strftime('%Y%m%d'), 2))
    folder = date.strftime("%d%b%y").lower()
    prefix = os.path.join(path, f'{folder}/SalishSea_1{res}_{datestr}')
    
    return prefix
colores=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def scatter_particles(ax, N ,n,nmin, nmax,yvar,lon,HD=0,colors='b'):
    '''scatter_particles(ax, N ,n,nmin, nmax,yvar,lon,HD=0,colors=colores)
    Use this function to scatter particles with different colours for each deploy location
    N= number of deploying sites,n=number of particles oper location, nmin,max=time min,max, yvar is the variable to plot on the yaxis, 
    Keep HD to 0 unless you want to plot with cartopy (only works for maps so yvar= latitude)
    '''
    
    scatter=[]
    #N is the number of stations, n the number of particles
    #Color is a list of strings picking the desired colors
    #nmin is t0, nmax is tmax
    #yvar is the y coordinate, lon is the longitud array
    #HD (1 for cartopy plot) 0 otherwise
    
    if HD == 0:
        if nmin==nmax:
            scatter.append(ax.scatter(lon[:, nmin], yvar[:, nmin],s=1,color=colors))
        else:
            
            scatter.append(ax.scatter(lon[:, nmin:nmax], yvar[:, nmin:nmax],s=1,color=colors))
    else:
        if nmin==nmax:
            scatter.append(ax.scatter(lon[:, nmin], yvar[:, nmin],s=1,transform=crs.PlateCarree(),zorder=2,color=colors))
        else:
            scatter.append(ax.scatter(lon[:, nmin:nmax], yvar[:, nmin:nmax],s=1,transform=crs.PlateCarree(),zorder=2,color=colors))
        
    return scatter

ss=[]

def mapanimation(outfile,N,n,clon,clat,fps=1,local=1):
    '''mapanimation(outfile,N,n,clon,clat,fps=1,local=1)
    Use this function to return an animated map of the particles,
    keep local=1 when working local and = 0 when remote. 
    outfile is the name of the output file from OP
    N= number of deploying sites,n=number of particles oper location,
    clat,clon location of deploying locations.
    '''
    coords,mask,ds = output(outfile,local)
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
    ss = []#scatter_particles(ax, N,n, 0,0, ds.lat,ds.lon)
    sed= {0: "w", 1: "k"}

    def update(frame):
        tstamp = ds.time[0, frame].values.astype('datetime64[s]').astype(datetime)
        t.set_text(tstamp.strftime('%Y-%b-%d %H:%M UTC'))
        global ss
        for scat in ss:
            scat.remove()
        ss = scatter_particles(ax, N,n, frame,frame, ds.lat,ds.lon)
        #ss.append(ax.scatter(ds.lon[:,frame], ds.lat[:,frame],c='m',s=5,alpha=ds.beached[:,frame].fillna(0))/3)
        #ss.append(ax.scatter(clon,clat,c='r', marker='*', linewidths=2))
        return ss
    return animation.FuncAnimation(fig, update, frames=np.arange(0,len(ds.lon[0,:]),fps))


def mapanimationd(outfile,N,n,clon,clat,fps=1,local=1):
    '''mapanimation(outfile,N,n,clon,clat,fps=1,local=1)
    Use this function to return an animated map of the particles,
    keep local=1 when working local and = 0 when remote. 
    outfile is the name of the output file from OP
    N= number of deploying sites,n=number of particles oper location,
    clat,clon location of deploying locations.
    '''
    coords,mask,ds = output(outfile,local)
    fig = plt.figure(figsize=(19, 8))
    ax = plt.axes(xlim=(-125,-122.5),ylim=(48.5,49.7))
    ax.contour(coords.nav_lon, coords.nav_lat, mask.mbathy[0,:,:],colors='k',linewidths=0.1)
    ax.contourf(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='lightgray')
    ax.contour(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='k')
    ax.set_aspect(1/1)
    dslo=len(ds.lon[0,:])
    plt.ylabel('Latitude',fontsize=16)
    plt.xlabel('Longitude',fontsize=16)
    t = ax.text(0.02, 0.02, '', transform=ax.transAxes)
    t.set_text('')
    cm = cmocean.cm.ice

    def update(frame):
        tstamp = ds.time[0, frame].values.astype('datetime64[s]').astype(datetime)
        t.set_text(tstamp.strftime('%Y-%b-%d %H:%M UTC'))
        
        ds2=ds.where(ds.time==ds.time[0,frame])
        ds2p=ds2.where(ds2.beached==0)
        dsb=ds2.where(ds2.beached==1)
        dss=ds2.where(ds2.beached==3)
        global ss
        for scat in ss:
            scat.remove()
        ss =[]
        
        
        #ss.append(ax.scatter(ds2.lon, ds2.lat,c='b',s=1))
        ss.append(ax.scatter(dsb.lon, dsb.lat,c='m',s=3))
        ss.append(ax.scatter(dss.lon, dss.lat,c='g',s=3))
        ss.append(ax.scatter(ds2p.lon, ds2p.lat,s=1,c=ds2p.tau,cmap=cm))
        
        print(f'{100*frame/dslo:.2f}% completed')
        return ss
    return animation.FuncAnimation(fig, update, frames=np.arange(0,len(ds.lon[0,:]),fps))


    

def visual(outfile,N,n,clon,clat,dmin,dd, nmin=0, nmax=-1,local=1):
    '''visual(outfile,N,n,clon,clat,dmin,dmax, nmin=0, nmax=-1,local=1)
    Use this function to return an animated map of the particles,
    keep local=1 when working local and = 0 when remote. 
    outfile is the name of the output file from OP
    N= number of deploying sites,n=number of particles oper location, dmin,dmax=deploying min,max depths,
    clat,clon location of deploying locations.
    '''
    coords,mask,ds = output(outfile,local)
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(19, 8))
    ax1.contour(coords.nav_lon, coords.nav_lat, mask.mbathy[0,:,:],colors='k',linewidths=0.1)
    ax1.contourf(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='lightgray')
    ax1.contour(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='k')

    scatter_particles(ax1, N,n, nmin, nmax, ds.lat,ds.lon)
    #ax1.scatter(clon,clat,c='g', marker='*', linewidths=1)

    scatter_particles(ax2, N,n, nmin, nmax, -ds.z,ds.lon)
    ax2.grid()
    plt.ylabel('Depth [m]')
    plt.xlabel('Longitude')

    if isinstance(dmin,int) :
        dmin=np.repeat((dmin),len(clon))
    else:
        dmin=[0-di for di in dmin]
    ax2.scatter(clon,-dmin,c='r', marker='*', linewidths=1)


def visuald(ax,outfile,N,n,clon,clat,dmin,dd, nmin=0, nmax=-1,local=1):
    '''visual(outfile,N,n,clon,clat,dmin,dmax, nmin=0, nmax=-1,local=1)
    Use this function to return an animated map of the particles,
    keep local=1 when working local and = 0 when remote. 
    outfile is the name of the output file from OP
    N= number of deploying sites,n=number of particles oper location, dmin,dmax=deploying min,max depths,
    clat,clon location of deploying locations.
    '''
    coords,mask,ds = output(outfile,local)
    #fig, (ax1, ax2) = plt.subplots(1,2,figsize=(19, 8))
    ax[0].contour(coords.nav_lon, coords.nav_lat, mask.mbathy[0,:,:],colors='k',linewidths=0.1)
    ax[0].contourf(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='lightgray')
    ax[0].contour(coords.nav_lon, coords.nav_lat, mask.tmask[0, 0, ...], levels=[-0.01, 0.01], colors='k')
   
    
    ds0=ds.where(ds.time<=ds.time[0,nmax])
    ds2=ds0.where(ds0.time>=ds0.time[0,nmin])
    dss=ds2.where(ds2.beached==3)
    scatter_particles(ax[0], N,n, 0, -1, ds2.lat,ds2.lon)
    ax[0].scatter(clon,clat,c='r', marker='*', linewidths=1)
    ax[0].set_ylabel('Longitude')
    ax[0].set_xlabel('Latitude')

    scatter_particles(ax[1], N,n, 0, -1, -ds2.z,ds2.lon)
    scatter_particles(ax[1], N,n, 0, -1, -dss.z,dss.lon,colors='r')
    t = ax[0].text(0.02, 0.02, '', transform=ax[0].transAxes)
    t.set_text('')
    tstamp = ds.time[0, nmax].values.astype('datetime64[s]').astype(datetime)
    t.set_text(tstamp.strftime('%Y-%b-%d %H:%M UTC'))
    ax[1].grid()
    ax[1].set_ylabel('Depth [m]')


    if isinstance(dmin,int) :
        dmin=np.repeat((dmin),len(clon))
    else:
        dmin=[0-di for di in dmin]
    ax[1].scatter(clon,-dmin,c='r', marker='*', linewidths=1)


def output(outfile,local=1):
    '''coords,mask,ds = output(outfile,local=1) 
    outfile is the name of the output file'''
    paths = path(local)
    coords = xr.open_dataset(paths['coords'], decode_times=False)
    mask = xr.open_dataset(paths['mask'])
    ds = xr.open_dataset(outfile)
    return  coords,mask,ds

labels0=['Nnm','Cmp','Vnc','Stl','Vct','Otf']
def profile(N,n,length,outfile,local=1,labels=labels0,levels=20,colors=colores):
    '''profile(N,n,length,outfile,levels=20,local=1)
    Use this function to return a depth profile of the particles,
    keep local=1 when working local and = 0 when remote. 
    outfile is the name of the output file from OP
    N= number of deploying sites,n=number of particles oper location, 
    length= number days of run,
    levels= how many layers to count particles
    '''
    coords,mask,ds = output(outfile,local)
    Z = np.linspace(0,430,levels)
    time = length*24+1
    zn = np.zeros([len(Z)-1,time])
    for j in range(time):
        zn[:,j],z_levels = np.histogram(ds.z[:, j], bins=Z)
    fig = plt.figure(figsize=(8, 8))
    ax = plt.axes(xlim=(-5,np.max(zn[:,0]+5)),ylim=(-500,0))
   
    plt.plot(zn[:,0],-z_levels[1:],'--')
    ax.grid()
    plt.ylabel('Depth [m]',fontsize=16)
    plt.xlabel('Particles',fontsize=16)

    plt.plot(zn[:,-1],-z_levels[1:],'-')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)

def profile2(ax,N,n,length,outfile,local=1,labels=labels0,levels=20,colors=colores):
    '''profile(N,n,length,outfile,levels=20,local=1)
    Use this function to return a depth profile of the particles,
    keep local=1 when working local and = 0 when remote. 
    outfile is the name of the output file from OP
    N= number of deploying sites,n=number of particles oper location, 
    length= number days of run,
    levels= how many layers to count particles
    '''
    coords,mask,ds = output(outfile,local)
    Z = np.linspace(0,430,levels)
    time = length*24+1
    zn = np.zeros([len(Z)-1,time])
    for j in range(time):
        zn[:,j],z_levels = np.histogram(ds.z[:, j], bins=Z)
   
    ax.plot(zn[:,0],-z_levels[1:],'--',label='$t_0$')
    ax.grid()
    ax.set_ylabel('Depth [m]')
    ax.set_title(outfile[30:39])
    ax.plot(zn[:,-1],-z_levels[1:],'-',label='$t_{end}$')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    ax.legend(fontsize=12)
    
    
def filename_set(start,length,varlist=['U','V','W'],local=0):
    '''filename,variables,dimensions = filename_set(start,duration,varlist=['U','V','W'],local=1)
    Modify function to include more default variables
    define start as: e.g, datetime(2018, 1, 17)
    length= number of days'''
    
    duration = timedelta(days=length)
    #Build filenames
    paths = path(local)
    Rlist,Tlist,Ulist, Vlist, Wlist = [], [], [], [], []
    Waveslist = []
    Flist = []
    Biolist,MZlist = [],[]
   
    for day in range(duration.days):
        path_NEMO = make_prefix(start + timedelta(days=day), paths['NEMO'])
        path_NEMO_d = make_prefix(start + timedelta(days=day), paths['NEMO'],res='d')
        Ulist.append(path_NEMO + '_grid_U.nc')
        Vlist.append(path_NEMO + '_grid_V.nc')
        Wlist.append(path_NEMO + '_grid_W.nc')
        Tlist.append(path_NEMO + '_grid_T.nc')
        Rlist.append(path_NEMO + '_carp_T.nc')
        Biolist.append(path_NEMO_d + '_prod_T.nc')
        MZlist.append(path_NEMO + '_ptrc_T.nc')
        Waveslist.append(get_WW3_path(start + timedelta(days=day)))
        Flist.append(get_Fraser_path(start + timedelta(days=day)))
        

    # Load NEMO forcing 
    filenames = {
        'U': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Ulist},
        'V': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Vlist},
        'W': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Wlist},
        'Kz': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Wlist},
        'T': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Tlist},
        'S': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Tlist},
        'R': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Rlist},
        'Bathy' : {'lon': paths['coords'], 'lat': paths['coords'], 'data': paths['mask']},
        'Cmask' : {'lon': paths['coords'], 'lat': paths['coords'],'depth': Wlist[0], 'data': paths['mask']},
        'D' : {'lon': paths['coords'], 'lat': paths['coords'], 'data': paths['mask']},
        'US' : {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist},
        'VS' : {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist},
        'WL' : {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist},
        'FS' :  {'lon': paths['coords'], 'lat': paths['coords'],'data': Flist},
        'MZ' : {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': MZlist},
        'Diat' : {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Biolist},
        'Flag' : {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Biolist},
    }
    variables = {'U': 'vozocrtx', 'V': 'vomecrty','W': 'vovecrtz','T':'votemper','S':'vosaline','R':'sigma_theta',
        'US':'uuss','VS':'vuss','WL':'lm','Bathy':'bathym', 'D':'Distc','FS':'rorunoff','Kz':'vert_eddy_diff',
        'MZ':'microzooplankton','Diat':'PPDIATNO3','Flag':'PPPHYNO3','Cmask':'coast_mask'}
    for fvar in varlist:
        if fvar == 'U' or fvar == 'Kz' or fvar == 'Diat' or fvar== 'Cmask':
            dimensions = {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw','time': 'time_counter'}
        elif fvar == 'US':
            dimensions = {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
        elif fvar == 'Bathy' or fvar == 'D' or fvar == 'FS':
            dimensions = {'lon': 'glamf', 'lat': 'gphif','time': 'time_counter'}
        
    
    file2,var2 = {},{}
    for var in varlist:
        file2[var]=filenames[var]
        var2[var]=variables[var]
    return file2,var2,dimensions

def p_deploy(N,n,dmin,dd,r = 1000):
    #r is radius of particle cloud [m]
    deg2m = 111000 * np.cos(50 * np.pi / 180)
    var = (r / (deg2m * 3))**2
    x_offset, y_offset = np.random.multivariate_normal([0, 0], [[var, 0], [0, var]], [n,N]).T
    if isinstance(dmin,int):
        zvals1 = dmin + np.random.random_sample([n,N]).T*(dd)
    else:
        zvals = []
        zvals1 = []
        for dept in dmin:
            zvals.append(dept + np.random.random_sample([n]).T*(dd))
        for i in range(len(zvals)):   
            zvals1=np.concatenate((zvals1[:],zvals[i]))
    return x_offset, y_offset, zvals1   


def get_WW3_path(date):
    """Construct WW3 results path given the date
    e.g., /opp/wwatch3/nowcast/SoG_ww3_fields_YYYYMMDD_YYYYMMDD.nc
    :arg date: date of WW3 record
    :type date: :py:class:`datetime.datetime`
    :returns: WW3 path
    :rtype: str
    """

    # Make WW3 path
    path = '/opp/wwatch3/nowcast'
    datestr = [date.strftime(fmt) for fmt in ('%d%b%y', '%Y%m%d_%Y%m%d')]
    path = os.path.join(path, datestr[0].lower(), f'SoG_ww3_fields_{datestr[1]}.nc')
    if not os.path.exists(path):
        raise ValueError(f"No WW3 record found for the specified date {date.strftime('%Y-%b-%d')}")

    return path


def p_unidist(lat0,lon0,bat,dy,dx):
    latbat = np.zeros(bat.shape)
    lonbat = np.zeros(bat.shape)

    for i in range(bat.shape[0]):
        for j in range(bat.shape[1]):
            if bat[i,j]>1e-5:
                latbat[i,j] = lat0[i,j]
                lonbat[i,j] = lon0[i,j]

    yi = np.arange(0,latbat.shape[0],dy)
    xi = np.arange(0,latbat.shape[1],dx)
    plat0,plon0 = latbat[yi,:],lonbat[yi,:]
    plat1,plon1 = plat0[:,xi],plon0[:,xi]

    plon,plat = [],[]
    for i in range(plat1.shape[0]):
        for j in range(plat1.shape[1]):
            if plat1[i,j]>1e-5:
                plat.append(plat1[i,j])
                plon.append(plon1[i,j])
    return plat,plon

def dist_coord(LAT,LON):
    la1,la2 = LAT[0],LAT[-1]
    lo1,lo2 = LON[0],LON[-1]
    R = 6378137
    PI=math.pi
    lat1 = la1 * PI/180 
    lat2 = la2 * PI/180
    dlat = (lat2-lat1)
    dlon = (lo1-lo2) * PI/180;
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = R * c
    return d/1e3


def p_unidist(lat0,lon0,bat,dy,dx):
    latbat = np.zeros(bat.shape)
    lonbat = np.zeros(bat.shape)

    for i in range(bat.shape[0]):
        for j in range(bat.shape[1]):
            if bat[i,j]>1e-5:
                latbat[i,j] = lat0[i,j]
                lonbat[i,j] = lon0[i,j]

    yi = np.arange(0,latbat.shape[0],dy)
    xi = np.arange(0,latbat.shape[1],dx)
    plat0,plon0 = latbat[yi,:],lonbat[yi,:]
    plat1,plon1 = plat0[:,xi],plon0[:,xi]

    plon,plat = [],[]
    for i in range(plat1.shape[0]):
        for j in range(plat1.shape[1]):
            if plat1[i,j]>1e-5:
                plat.append(plat1[i,j])
                plon.append(plon1[i,j])
    return plat,plon

#def load_config(config_yaml):
#   with open(config_yaml) as f:
#       config = yaml.safe_load(f)
#   return config


def load_config(config_yaml):
   with open(config_yaml[0]) as f:
       config = yaml.safe_load(f)
   return config


def get_Fraser_path(date):
    """Construct Fraser river outflow path given the date
    e.g., /results/forcing/rivers/#R201702DFraCElse_yYYYYmMMdDD.nc
    :arg date: date of fraser outflow record 
    :type date: :py:class:`datetime.datetime`
    :returns: nc path
    :rtype: str
    """

    # Make Fraser path
    path = '/results/forcing/rivers/'
    datestr = [date.strftime(fmt) for fmt in ('%d%b%y', 'y%Ym%md%d')]
    path = os.path.join(path, f'R201702DFraCElse_{datestr[1]}.nc')
    if not os.path.exists(path):
        raise ValueError(f"No file found for the specified date {date.strftime('%Y-%b-%d')}")
    return path

def get_timestamps(start,length):
    timestamps=[]
    duration = timedelta(days=length)
    for day in range(duration.days):
        timestamps.append([start + timedelta(days=day)])
    return np.array(timestamps, dtype='datetime64')

def kernel_asem(pset,config):
    KER = AdvectionRK4_3D 
    if 'Buoyancy' in config['kernel']:
        KER += pset.Kernel(Buoyancy)
    if 'Stokes_drift' in config['kernel']:
        KER += pset.Kernel(Stokes_drift)
    if 'Beaching' in config['kernel']:
        KER += pset.Kernel(Beaching)
        KER += pset.Kernel(Unbeaching)
    if 'Turb_mix' in config['kernel']:
        KER += pset.Kernel(turb_mix)
    if 'Turb_mix2' in config['kernel']:
        KER += pset.Kernel(turb_mix2)
    if 'Biofilm' in config['kernel']:
        KER += pset.Kernel(Biofilm)
    return KER

def particle_maker(config):
    #Define particle properties 
    class MPParticle(JITParticle):
        if 'ro' in config['particle']:        
            ro = Variable('ro', initial = config['particle']['ro'])  # config['particle']['ro']
        if 'diameter' in config['particle']:           
            diameter = Variable('diameter', initial = config['particle']['diameter'])
        if 'SDD' in config['particle']:  
            Sdd = Variable('SDD', initial = config['particle']['SDD'])
        if 'SDL' in config['particle']:           
            Sdl = Variable('SDL', initial = config['particle']['SDL'])
        if 'length' in config['particle']:  
            length = Variable('length', initial = config['particle']['length'])
        if 'Lb' in config['particle']:  
            Lb = Variable('Lb', initial = config['particle']['Lb'])  #days needed in days for particle to have 67% probability of beaching if in beaching zone (500m)
        if 'Db' in config['particle']:  
            Db = Variable('Db', initial = config['particle']['Db']) #Distance at which particles can randomly beach.
        if 'Ub' in config['particle']:  
            Ub = Variable('Ub', initial = config['particle']['Ub'])  #days to have 67% probability of unbeaching
        if 'beached' in config['particle']:  
            beached = Variable('beached', initial = 0)
        if 'Ws' in config['particle']:  
            Ws = Variable('Ws', initial =  config['particle']['Ws']) #200m/dia
        if 'tau' in config['particle']:  
            tau = Variable('tau', initial =  config['particle']['tau']) # track number of particles
        if 'fratio' in config['particle']:  
            fratio = Variable('fratio', initial =  config['particle']['fratio']) # track number of particles
        if 'Nbac' in config['particle']:  
            Nbac = Variable('Nbac', initial =  0) # number of bacteria attached
        if 'Nflag' in config['particle']:  
            Nflag = Variable('Nflag', initial =  0) # number of flagellates grazing on the attached bacteria
        if 'dz' in config['particle']:  
            dz = Variable('dz', initial =  0) # dz variable
    return MPParticle
