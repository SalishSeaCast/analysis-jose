import sys 
import xarray as xr
import numpy as np
import os
import pandas as pd
import yaml
import math
from datetime import datetime, timedelta
from parcels import FieldSet, Field, VectorField, ParticleSet, JITParticle, ParcelsRandom, Variable,StatusCode,GeographicPolar, Geographic
from glob import glob
import warnings
from mpi4py import MPI
import time, resource, json
from types import SimpleNamespace
import shutil

comm = MPI.COMM_WORLD; rank = comm.Get_rank()

warnings.filterwarnings("ignore")

sys.path.append('/home/jvalenti/MOAD/analysis-jose/notebooks/MP_model/Source')

from OP_Kernels import *
from OP_functions import *


from scipy.spatial import cKDTree



def _morton(ty, tx, bits=16):
    ty = ty.astype(np.uint64); tx = tx.astype(np.uint64)
    code = np.zeros(len(ty), dtype=np.uint64)
    one = np.uint64(1)
    for b in range(bits):
        sh = np.uint64(b)
        code |= ((tx >> sh) & one) << np.uint64(2 * b)
        code |= ((ty >> sh) & one) << np.uint64(2 * b + 1)
    return code


def make_chunk_aware_partition(lon_grid, lat_grid, chunksize=None, tile=None,
                               coords_order='lonlat'):
    """Equal-count, spatially compact partition aligned to dask chunks.
    tile: tile size in grid cells (overrides chunksize, needed for 'auto').
    coords_order: 'lonlat' or 'latlon', the column order Parcels passes in coords."""
    ny, nx = lon_grid.shape
    if tile is not None:
        tile_y = tile_x = int(tile)
    else:
        tile_y, tile_x = int(chunksize['lat'][1]), int(chunksize['lon'][1])

    valid = np.isfinite(lon_grid) & np.isfinite(lat_grid)
    flat_idx = np.flatnonzero(valid)
    cos_lat = np.cos(np.deg2rad(np.nanmean(lat_grid[valid])))
    tree = cKDTree(np.column_stack((lon_grid[valid] * cos_lat, lat_grid[valid])))

    def partition_function(coords, mpi_size=1):
        n = len(coords)
        if mpi_size == 1:
            return np.zeros(n, dtype=np.int32)
        if coords_order == 'lonlat':
            plon, plat = coords[:, 0], coords[:, 1]
        else:
            plat, plon = coords[:, 0], coords[:, 1]

        _, nearest = tree.query(np.column_stack((plon * cos_lat, plat)))
        cell = flat_idx[nearest]
        j, i = np.unravel_index(cell, (ny, nx))
        key = _morton(j // tile_y, i // tile_x)

        # primary: tile (Morton), then cell, then particle index -> deterministic
        order = np.lexsort((np.arange(n), cell, key))
        blocks = (np.arange(n, dtype=np.int64) * mpi_size) // n   # 0..mpi_size-1, sizes differ by <= 1
        out = np.empty(n, dtype=np.int32)
        out[order] = blocks
        return out

    return partition_function


def add_constants(field_set,param):
    field_set.add_constant('deg2met', param['constants']['deg2met'])
    field_set.add_constant('deglatT', param['constants']['deglatT'])   
    field_set.add_constant('kh', param['constants']['kh'])  
    field_set.add_constant('rho_hb', param['constants']['rho_hb']) 
    field_set.add_constant('rho_ph', param['constants']['rho_ph']) 
    field_set.add_constant('V_hb', param['constants']['V_hb'])
    field_set.add_constant('V_ph', param['constants']['V_ph']) 
    field_set.add_constant('gravity', param['constants']['gravity']) 
    field_set.add_constant('biol_dt', param['constants']['biol_dt']) 
    field_set.add_constant('dyn_dt', param['constants']['dyn_dt']) 
    field_set.add_constant('kB', param['constants']['kB']) 
    field_set.add_constant('N_hb_cell', param['constants']['N_hb_cell'])
    field_set.add_constant('N_ph_cell', param['constants']['N_ph_cell']) 
    field_set.add_constant('r_hbac', param['constants']['r_hbac']) 
    field_set.add_constant('r_phyt', param['constants']['r_phyt']) 
    field_set.add_constant('vs_hbac', param['constants']['vs_hbac']) 
    field_set.add_constant('vs_diat', param['constants']['vs_diat']) 
    field_set.add_constant('resusp_h', param['constants']['resusp_h']) 
    PrUb = 1- np.exp(-param['param']['dt']/(param['constants']['Ub']*86400))
    field_set.add_constant('PrUb', PrUb)

    return field_set

def build_fieldset(start,duration,cs):
    paths = path()
    cgrid_dim = {'lon': 'glamf','lat': 'gphif','depth':'depthw', 'time':'time_counter'}
    Ulist, Vlist, Wlist, Tlist, Blist, Plist, Glist ,Waveslist, Windlist = [], [], [], [], [], [], [], [], []

    for day in range(duration.days):
        path_NEMO = make_prefix(start + timedelta(days=day),paths['NEMO'])
        Ulist.append(path_NEMO + '_grid_U.nc')
        Vlist.append(path_NEMO + '_grid_V.nc')
        Wlist.append(path_NEMO + '_grid_W.nc')
        Tlist.append(path_NEMO + '_grid_T.nc')

        Glist.append(get_Biol_path(start + timedelta(days=day),'graz'))
        Plist.append(get_Biol_path(start + timedelta(days=day),'prod'))
        Blist.append(get_Biol_path(start + timedelta(days=day)))
        
        Waveslist.append(get_WW3_path(start + timedelta(days=day)))
        Windlist.append(get_Winds_path(start + timedelta(days=day)))

    filenames = {
    'U': {'lon': paths['coords'],'lat': paths['coords'],'depth': Wlist[0],'data': Ulist},
    'V': {'lon': paths['coords'],'lat': paths['coords'],'depth': Wlist[0],'data': Vlist},
    'W': {'lon': paths['coords'],'lat': paths['coords'],'depth': Wlist[0],'data': Wlist}
    }
    dimensions = {
    'U': {'lon': 'glamf','lat': 'gphif','depth':'depthw', 'time':'time_counter'},
    'V': {'lon': 'glamf','lat': 'gphif','depth':'depthw', 'time':'time_counter'},
    'W': {'lon': 'glamf','lat': 'gphif','depth':'depthw', 'time':'time_counter'},
    }
    variables = {
    'U': 'vozocrtx',
    'V': 'vomecrty',
    'W': 'vovecrtz'
    }

    field_set = FieldSet.from_nemo(
        filenames,
        variables,
        dimensions,
        allow_time_extrapolation=True,
        chunksize=cs
    )

    Kz = Field.from_netcdf(filenames['W'], 'vert_eddy_diff', dimensions['W'],allow_time_extrapolation=True, chunksize=cs)
    field_set.add_field(Kz)

    filenames['tracer'] = {'lon': paths['coords'],'lat': paths['coords'],'depth': Tlist[0],'data': Tlist}
    dimensions['tracer'] = {'lon': 'glamt', 'lat': 'gphit', 'depth': 'deptht','time': 'time_counter'}

    R = Field.from_netcdf(filenames['tracer'], 'sigma_theta', dimensions['tracer'],allow_time_extrapolation=True, chunksize=cs)
    T = Field.from_netcdf(filenames['tracer'], 'votemper', dimensions['tracer'],allow_time_extrapolation=True, chunksize=cs)

    field_set.add_field(R)
    field_set.add_field(T)

    #BIOLOGY fields
    filenames['Biology'] = {'lon': paths['coords'],'lat': paths['coords'],'depth': Blist[0],'data': Blist}
    filenames['Prod'] = {'lon': paths['coords'],'lat': paths['coords'],'depth': Plist[0],'data': Plist}
    filenames['Mort'] = {'lon': paths['coords'],'lat': paths['coords'],'depth': Glist[0],'data': Glist}

    PP = Field.from_netcdf(filenames['Prod'], 'PPDIAT', dimensions['tracer'],allow_time_extrapolation=True, chunksize=cs)
    HBProd = Field.from_netcdf(filenames['Prod'], 'HetHBAC', dimensions['tracer'],allow_time_extrapolation=True, chunksize=cs)
    Diat_mort = Field.from_netcdf(filenames['Mort'], 'MORTDIAT', dimensions['tracer'],allow_time_extrapolation=True, chunksize=cs)
    Hbac_mort = Field.from_netcdf(filenames['Mort'], 'MORTHBAC', dimensions['tracer'],allow_time_extrapolation=True, chunksize=cs)
    Hbac = Field.from_netcdf(filenames['Biology'], 'heterotrophic_bacteria', dimensions['tracer'],allow_time_extrapolation=True, chunksize=cs)
    Diat = Field.from_netcdf(filenames['Biology'], 'diatoms', dimensions['tracer'],allow_time_extrapolation=True, chunksize=cs)

    field_set.add_field(PP)
    field_set.add_field(HBProd)
    field_set.add_field(Diat_mort)
    field_set.add_field(Hbac_mort)
    field_set.add_field(Hbac)
    field_set.add_field(Diat)

    #####FROM HERE ALL NOT CHUNKED FIELDS (2D or different dimensions)
    #Add Stokes Drift fields
    dimensions['stokes'] = {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
    filenames['stokes'] = {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist}
    us = Field.from_netcdf(filenames['stokes'],'uuss', dimensions['stokes'],allow_time_extrapolation=True)
    vs = Field.from_netcdf(filenames['stokes'], 'vuss', dimensions['stokes'],allow_time_extrapolation=True)
    wl = Field.from_netcdf(filenames['stokes'], 'lm', dimensions['stokes'],allow_time_extrapolation=True)
    us.units = GeographicPolar()   # u: m/s -> deg lon/s 
    vs.units = Geographic() 
    field_set.add_field(us)
    field_set.add_field(vs)
    field_set.add_field(wl)
    field_set.add_vector_field(VectorField("stokes", us, vs, wl))

    filenames['totaldepth'] = {'lon': paths['coords'], 'lat': paths['coords'], 'data': paths['mask']}
    filenames['ssh'] = {'lon': paths['coords'], 'lat': paths['coords'], 'data': Tlist}

    #Add Bathymetry 2D field
    dimensions['flat'] = {'lon': 'glamt', 'lat': 'gphit'}
    TD = Field.from_netcdf(filenames['totaldepth'], 'totaldepth', dimensions['flat'],allow_time_extrapolation=True, interp_method='nearest')
    #TD = Field.from_netcdf(filenames['totaldepth'], 'totaldepth', dimensions['flat'],allow_time_extrapolation=True)
    field_set.add_field(TD)

    #Add SSH 
    dimensions['flat_t'] = {'lon': 'glamt', 'lat': 'gphit','time': 'time_counter'}
    SSH = Field.from_netcdf(filenames['ssh'], 'sossheig', dimensions['flat_t'],allow_time_extrapolation=True)
    field_set.add_field(SSH)

    b_stress = Field.from_netcdf(filenames['Biology'], 'bottom_stress', dimensions['flat_t'],allow_time_extrapolation=True)
    field_set.add_field(b_stress)

    filenames['wind'] = {'lon': Windlist[0],'lat': Windlist[0],'data': Windlist}
    dimensions['wind'] = {'lon': 'nav_lon', 'lat': 'nav_lat','time': 'time_counter'}

    u_wind = Field.from_netcdf(filenames['wind'], 'u_wind', dimensions['wind'],allow_time_extrapolation=True)
    v_wind = Field.from_netcdf(filenames['wind'], 'v_wind', dimensions['wind'],allow_time_extrapolation=True)
    u_wind.units = GeographicPolar()   # u: m/s -> deg lon/s 
    v_wind.units = Geographic()  
    field_set.add_field(u_wind)
    field_set.add_field(v_wind)
    field_set.add_vector_field(VectorField("wind", u_wind, v_wind))
    return field_set


    
def mpi_Parcels(config='config.yaml', runtime=None, outputdt=None,
             outfile=None, particles=None, extra_kernels=(), convert=True, seed=None, cs_override=None, overwrite=False):
    
    param = load_config(config)
    start = datetime(param['startdate']['year'], param['startdate']['month'], param['startdate']['day'])
    duration = timedelta(days=param['param']['length'])   # sets which files get loaded

    runtime_h = runtime.total_seconds()/3600 if runtime is not None else duration.days*24

    # particles: None -> read from deploy file; dict -> use as given
    if particles is None:
        MP = particle_source(param, seed=param['param'].get('seed', 1234), runtime_hours=runtime_h)
    else:
        MP = SimpleNamespace(**particles)
    lat, lon, z = MP.lat, MP.lon, MP.depth
    a, b, c, Vol_0, density, SI, Time = MP.a, MP.b, MP.c, MP.Vol, MP.rho_p, MP.SI, MP.time
    
    cs = {
        "time": ("time", param['cs']['time']),
        "depth": ("depth", param['cs']['depth']),   
        "lat": ("latitude", param['cs']['lat']),
        "lon": ("longitude", param['cs']['lon']),
    } #fieldsetchunking
    
    if cs_override == 'none':
        cs = None
    elif cs_override is not None:
        cs = cs_override
    cs_tag = ('none' if cs is None else 'auto' if isinstance(cs, str)
            else f"y{cs['lat'][1]}_x{cs['lon'][1]}_z{cs['depth'][1]}")

    dt = param['param']['dt'] #mixing timestep 
    dtp = param['param']['dtp'] 
    repeatdt = timedelta(hours=dtp) if dtp > 0 else None
    odt = param['param']['odt']

    if rank == 0:
        print(f"The model will run for {duration.days} days, starting at {start}")

    if outfile is None:
        outfile = os.path.join('/home/jvalenti/MOAD/results/', param['file']['name'] + '.zarr')
    static_file = outfile.replace('.zarr', '_static.nc')

    if rank == 0:
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
        for p in (outfile, static_file):
            if os.path.exists(p):
                if not overwrite:
                    raise FileExistsError(f"{p} exists; pass overwrite=True to replace it")
                shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
    comm.Barrier()

    t0 = time.perf_counter()
    field_set = build_fieldset(start,duration,cs)
    t_build = time.perf_counter() - t0

    field_set = add_constants(field_set,param)
    #Define our Parcels Particles
    class MPParticle(JITParticle):    
        fact = Variable('fact', initial = 1.0)
        status = Variable('status', initial = 1) #0=non initialized, 1=active, 2=beached, 3 = sedimented, 5=export to PO, 10=lost
        N_hb = Variable('N_hb', initial = 0)  #n of hbac attached
        N_ph = Variable('N_ph', initial = 0)  #n of diat attached
        bioth = Variable('bioth', initial = 0.0, to_write=False)  #thickness biofouling (can be calculated from N and mp dimensions)
        #these do not change change in run
        dim1 = Variable('dim1', initial = param['particle']['dim1'], to_write=False) # particle small semiaxis
        dim2= Variable('dim2', initial = param['particle']['dim2'], to_write=False)  # particle medium semiaxis
        dim3 = Variable('dim3', initial = param['particle']['dim3'], to_write=False) # particle large semiaxis
        density = Variable('density', initial = param['particle']['density'], to_write=False) #particle density 
        shape = Variable('shape', initial = param['particle']['shape'], to_write=False) # particle roundness (0 is cube, 1 is sphere) 
        Vol = Variable('Vol', initial = 1e-13, to_write=False)  #estimated at the initialization
        ws = Variable('ws', initial=0.0)
        last_update_biol = Variable('last_update_biol', initial=-5000.0, to_write=False) #careful with initialization
        last_update_dyn = Variable('last_update_dyn', initial=-100.0, to_write=False)
        vtau_crit = Variable('vtau_crit', initial=1.0, to_write=False)
    
    static = xr.Dataset(
    {name: ('release_index', np.asarray(getattr(MP, name), dtype=np.float64))
     for name in ['a', 'b', 'c', 'Vol', 'rho_p', 'SI', 'lon', 'lat', 'depth']}, coords={'release_index': np.arange(len(MP.lat))},
    attrs={'N_per_release': len(MP.lat), 'repeatdt_hours': dtp})

    if rank == 0:
        static.to_netcdf(static_file)
    comm.Barrier()
   
    if seed is not None:
        ParcelsRandom.seed(seed + rank)

    co = xr.open_dataset(path()['coords'], decode_times=False)
    lon_t, lat_t = co.glamt.squeeze().values, co.gphit.squeeze().values

    partition = make_chunk_aware_partition(
        lon_grid=lon_t, lat_grid=lat_t,
        chunksize=cs,tile=200)          # see point 3 

    pset = ParticleSet.from_list(field_set, MPParticle, lon=lon, lat=lat, depth=z, dim3 = a, dim2 = b, dim1 = c, Vol = Vol_0, density= density,time = Time,  partition_function=partition)

    t0 = time.perf_counter()
    
    pset.execute(
        [vert_mix,Advection,horz_mix,Stokes_Windage,LandCheck,Biofouling_ws,Buoyancy,export,KeepInOcean,CheckOutOfBounds,*extra_kernels],
        runtime=runtime if runtime is not None else duration,
        dt=timedelta(seconds=dt),
        output_file=pset.ParticleFile(name=outfile, outputdt=timedelta(hours=odt),chunks=(int(1e4), 1)))

    t_exec = time.perf_counter() - t0

    mem_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024   # KB -> MB on Linux
    t_all = comm.gather(t_exec, root=0)
    m_all = comm.gather(mem_mb, root=0)
    comm.Barrier()

    n_local = len(pset)
    t_all = comm.gather(t_exec, root=0)
    m_all = comm.gather(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, root=0)
    n_all = comm.gather(n_local, root=0)
    if rank == 0:
        print("BENCH " + json.dumps(dict(
            ranks=comm.Get_size(), N=len(lat), cs=cs_tag, t_build=t_build,
            t_exec_max=max(t_all), t_exec_min=min(t_all),
            n_per_rank_min=min(n_all), n_per_rank_max=max(n_all),
            mem_max_mb=max(m_all))))

    if convert and rank == 0:
        ds, name = zarr_tonet(outfile)
        ds.to_netcdf('/home/jvalenti/scratch/Parcels/results/' + name + '.nc')

    

if __name__=="__main__":
    try:
        config = sys.argv[1]
    except :
        print('Something went wrong')
    mpi_Parcels(config)
