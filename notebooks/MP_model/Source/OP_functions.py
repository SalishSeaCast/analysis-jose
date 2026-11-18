import pandas as pd
import numpy as np
import yaml 
from datetime import datetime, timedelta
import xarray as xr
import os
from glob import glob
import math
from types import SimpleNamespace

def path():
    machine = os.uname().nodename
    F = os.environ.get('FORCING_DIR')   
    if 'salish' in machine:
        paths = {'NEMO': '/results2/SalishSea/nowcast-green.202111',
            'coords': '/ocean/jvalenti/MOAD/grid/coordinates_seagrid_SalishSea201702.nc',
            'coordsWW3': '/ocean/jvalenti/MOAD/grid2/WW3_grid.nc',
            'mask': '/ocean/jvalenti/MOAD/grid2/mesh_mask202108_TDV.nc',
            'bat': '/ocean/jvalenti/MOAD/grid/bathymetry_202108.nc',
            'out': '/home/jvalenti/MOAD/results',
            'WW3': '/opp/wwatch3/',
            'SHEM': '/data/jvalenti/results',
            'winds': '/results/forcing/atmospheric/continental2.5/nemo_forcing'}
    elif 'nibi' in machine:
        paths = {'NEMO': '/project/def-allen/SalishSea/nowcast-green.202111/',
            'coords': '/home/jvalenti/MOAD/grid/coordinates_seagrid_SalishSea201702.nc',
            'coordsWW3': '/home/jvalenti/MOAD/grid2/WW3_grid.nc',
            'mask': '/home/jvalenti/MOAD/grid2/mesh_mask202108_TDV.nc',
            'bat': '/home/jvalenti/MOAD/grid/bathymetry_202108.nc',
            'out': '/scratch/jvalenti/MP_runs',
            'home': '/home/jvalenti/MOAD/analysis-jose/'}
    elif any(name in machine for name in ['narval', 'nc']):
        paths = {'NEMO': '/home/jvalenti/scratch/Parcels/forcing/SalishSeaCast',
            'coords': '/home/jvalenti/MOAD/grid/coordinates_seagrid_SalishSea201702.nc',
            'coordsWW3': '/home/jvalenti/MOAD/grid2/WW3_grid.nc',
            'mask': '/home/jvalenti/MOAD/grid2/mesh_mask202108_TDV.nc',
            'bat': '/home/jvalenti/MOAD/grid/bathymetry_202108.nc',
            'out': '/home/jvalenti/scratch/Parcels/results',
            'WW3': '/home/jvalenti/scratch/Parcels/forcing/WW3',
            'SHEM': '/home/jvalenti/scratch/Parcels/forcing/SHEM',
            'winds': '/home/jvalenti/projects/def-allen/SalishSea/forcing/forcing/atmospheric/continental2.5/nemo_forcing'}
    if F:
            paths['NEMO'] = os.path.join(F, 'SalishSeaCast')
            paths['SHEM'] = os.path.join(F, 'SHEM')
            paths['WW3'] = os.path.join(F, 'WW3')
    return paths

def get_Biol_path(date,filen='biol'):
    # Make biol path
    pathshem = path()['SHEM']
    datestr = [date.strftime(fmt) for fmt in ('%d%b%y', '%Y%m%d-%Y%m%d')]
    if filen =='biol':
        pathshem = os.path.join(pathshem, datestr[0].lower(), f'SHEM_1h_{datestr[1]}_biol_T.nc')
    elif filen =='prod':
        pathshem = os.path.join(pathshem, datestr[0].lower(), f'SHEM_1d_{datestr[1]}_prod_T.nc')
    elif filen =='graz':
            pathshem = os.path.join(pathshem, datestr[0].lower(), f'SHEM_1d_{datestr[1]}_graz_T.nc')
    return pathshem  

def get_Winds_path(date):
    """Construct Winds results path given the date
    :rtype: str
    """
    datestr = date.strftime('y%Ym%md%d')
    # Make winds path
    if date > datetime(2023,2,22):
        pathwinds = path()['winds']
        pathwinds = os.path.join(pathwinds, f'hrdps_{datestr}.nc')
    else:
        pathwinds = '/results/forcing/atmospheric/GEM2.5/operational'
        pathwinds = os.path.join(pathwinds, f'ops_{datestr}.nc')
    return pathwinds  

def get_WW3_path(date):
    """Construct WW3 results path given the date
    e.g., /opp/wwatch3/nowcast/SoG_ww3_fields_YYYYMMDD_YYYYMMDD.nc
    :arg date: date of WW3 record
    :type date: :py:class:`datetime.datetime`
    :returns: WW3 path
    :rtype: str
    """
    # Make WW3 path
    machine = os.uname().nodename
    if 'salish' in machine:
        path0 = path()['WW3']
        path1 = path0+ 'hindcast'
        path2 = path0+ 'nowcast'
        datestr = [date.strftime(fmt) for fmt in ('%d%b%y', '%Y%m%d_%Y%m%d')]
        path1 = os.path.join(path1, datestr[0].lower(), f'SoG_ww3_fields_{datestr[1]}.nc')
        if not os.path.exists(path1):
            path2 = os.path.join(path2, datestr[0].lower(), f'SoG_ww3_fields_{datestr[1]}.nc')
            if not os.path.exists(path2):    
                raise ValueError(f"No WW3 record found for the specified date {date.strftime('%Y-%b-%d')}")
    elif any(name in machine for name in ['narval', 'nc']):
        path0 = path()['WW3']
        datestr = [date.strftime(fmt) for fmt in ('%d%b%y', '%Y%m%d_%Y%m%d')]
        path0 = os.path.join(path0, datestr[0].lower(), f'SoG_ww3_fields_{datestr[1]}.nc')
        if not os.path.exists(path0):    
            raise ValueError(f"No WW3 record found for the specified date {date.strftime('%Y-%b-%d')}")
    return path0  

def load_config(config_yaml):
   with open(config_yaml) as f:
       config = yaml.safe_load(f)
   return config

def zarr_tonet(fileoutname):
    name = os.path.splitext(os.path.basename(fileoutname))[0]
    files = sorted(glob(os.path.join(fileoutname, "proc*")))
    if files:   # MPI run: one store per rank
        ds = xr.concat(
            [xr.open_zarr(f) for f in files],
            dim="trajectory",
            compat="no_conflicts",
            coords="minimal",
        )
    else:       # serial run: single store
        ds = xr.open_zarr(fileoutname)
    return ds, name

def simple_partition_function(coords, mpi_size=1):
    """A very simple partition function
    that assigns particles to processors
    """
    return np.linspace(0, mpi_size, coords.shape[0], endpoint=False, dtype=np.int32)


def make_prefix(date, path, res='h'):
    """Construct path prefix for local SalishSeaCast results given date object and paths dict
    e.g., /results2/SalishSea/nowcast-green.201905/daymonthyear/SalishSea_1h_yyyymmdd_yyyymmdd
    """
    datestr = '_'.join(np.repeat(date.strftime('%Y%m%d'), 2))
    folder = date.strftime("%d%b%y").lower()
    prefix = os.path.join(path, f'{folder}/SalishSea_1{res}_{datestr}')
    return prefix

def pandas_deploy(N,MFc,r,dd,dtp):
    n =1
    MFc = float(MFc)
    Outfall_deploy = pd.read_csv(N, index_col = [0])
    Pol = list(Outfall_deploy.Population)
    Lat = Outfall_deploy.Lat
    Lon = Outfall_deploy.Lon
    Depth = Outfall_deploy.Depth
    clat = []
    clon = []
    cz = []
    for i,loc in enumerate(Pol):
        for j in range(int(round((loc*250*dtp)/MFc,0))):
            clat.append(Lat.iat[i])
            clon.append(Lon.iat[i])
            try:
                cz.append(float(Depth.iat[i]))
            except ValueError: #This is for rivers
                cz.append(2)
    N = len(clat)
    deg2m = 111319.5 * np.cos(49 * np.pi / 180)
    var = (r / (deg2m * 3))**2
    x_offset, y_offset = np.random.multivariate_normal([0, 0], [[var, 0], [0, var]], [n,N]).T #Here we define a cloud of deployment points.
    z_offset=np.random.random_sample([N]).T*(dd) #add some depth variability
    lon = np.zeros([N])
    lat = np.zeros([N])
    z = np.zeros([N])
    for i in range(N):
        lon[i]=(clon[i] + x_offset[i])
        lat[i]=(clat[i] + y_offset[i])
        z[i]=(cz[i] + z_offset[i])
    return lon,lat,z

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


def particle_source(param, seed=None,runtime_hours=None):
    """Build particle release arrays from the deployment file

    Each row releases round(number / SI) super-individuals, where
    SI = param['constants']['SI']. Rows that round to 0 are skipped.

    Returns a namespace with float64 arrays (one entry per particle):
        lat, lon, depth, a, b, c, Vol, rho_p
    with a >= b >= c (a = large semiaxis, c = small), matching the
    dim3/dim2/dim1 convention in MPParticle and the kernels.
    """
    rng = np.random.default_rng(seed)
    P = param['particle']
    C = param['constants']
    dtp = param['param']['dtp']                      # hours between releases; 0 = single pulse
    if runtime_hours is None:
        runtime_hours = param['param']['length'] * 24.0

    # ---- read deployment file ----
    df = pd.read_csv(P['filename'])
    lat0 = df.iloc[:, 2].to_numpy(dtype=float)
    lon0 = df.iloc[:, 3].to_numpy(dtype=float)
    z0 = pd.to_numeric(df.iloc[:, 5], errors='coerce').fillna(2.0).to_numpy(dtype=float)
    npart = df.iloc[:, 4].to_numpy(dtype=float)*250*runtime_hours #particles in run

    # super-individuals per location
    # weight npart % SI if there is a remainder.
    SI0 = int(C['SI'])
    npart = npart.astype(np.int64)      
    n_full, rem = np.divmod(npart, SI0)
    counts = n_full + (rem > 0)   
    idx = np.repeat(np.arange(len(npart)), counts)  
    SI = np.full(len(idx), float(SI0))
    last = np.cumsum(counts) - 1
    SI[last[rem > 0]] = rem[rem > 0]  

    clat, clon, cz = lat0[idx], lon0[idx], z0[idx]
    N = len(clat)
    if N == 0:
        raise ValueError(f"No particles to release: check {P['filename']}")

     # release time (s): stratified over the run within each location
    j = np.arange(N) - (np.cumsum(counts) - counts)[idx]      # 0..counts-1 within row
    u = rng.uniform(0.0, 1.0, N)
    u[j == 0] = 0.0                                           # first SI of each row at t = 0
    t_h = (j + u) / counts[idx] * runtime_hours
    if dtp > 0:
        t_h = np.floor(t_h / dtp) * dtp                       # snap to the release grid
    else:
        t_h = np.zeros(N)                                     # single pulse
    time = t_h * 3600.0

    # cloud of release points
    sig_lat = (P['r'] / 3.0) / C['deg2met']    
    sig_lon = (P['r'] / 3.0) / C['deglatT']      
    lon = clon + rng.normal(0.0, sig_lon, N)
    lat = clat + rng.normal(0.0, sig_lat, N)
    z = cz + rng.uniform(0.0, P['dd'], N)  # extra depth in [0, dd]

    # particle dimensions (semiaxes, m) floor at 10% of the mean,
    def dimension_sd(mean, sd):
        return np.clip(rng.normal(mean, sd, N), 0.1 * mean, None)

    c = dimension_sd(P['dim1'], P.get('dim1sd', 0.0))
    b = dimension_sd(P['dim2'], P.get('dim2sd', 0.0))
    a = dimension_sd(P['dim3'], P.get('dim3sd', 0.0))

    # superellipsoid volume
    e = P['shape']
    fe = e**2 * math.gamma(e/2)**3 / math.gamma(3*e/2) 
    Vol = (2.0/3.0)*fe*a*b*c

    rho_p = np.full(N, P['density'], dtype=float)

    print(f"{N} super-individuals over {runtime_hours:g} h "
              f"({len(lat0)} locations, SI = {C['SI']:g}, sum SI = {SI.sum():.0f})")

    return SimpleNamespace(lat=lat, lon=lon, depth=z, a=a, b=b, c=c,
                           Vol=Vol, rho_p=rho_p, SI=SI, time=time)

