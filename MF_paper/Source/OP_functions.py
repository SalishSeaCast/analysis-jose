def load_config(config_yaml):
   with open(config_yaml[0]) as f:
       config = yaml.safe_load(f)
   return config

def zarr_tonet(fileoutname):
    from os import path
    name =  fileoutname.split('_')[1]
    fname = fileoutname
    print(fileoutname)
    files = glob(path.join(fname, "proc*"))
    ds = xr.concat(
        [xr.open_zarr(f) for f in files],
        dim="trajectory",
        compat="no_conflicts",
        coords="minimal",
    )
    return ds,name

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
    prefix = f'{path}/SalishSea_1{res}_{datestr}'
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

def filename_set(start,length,varlist=['U','V','W']):
    '''filename,variables,dimensions = filename_set(start,duration,varlist=['U','V','W'],local=1)
    Modify function to include more default variables
    define start as: e.g, datetime(2018, 1, 17)
    length= number of days'''
    
    duration = timedelta(days=length)
    #Build filenames
    paths = path()
    Tlist,Ulist, Vlist, Wlist = [], [], [], [], []
    Waveslist = []
   
    for day in range(duration.days):
        path_NEMO = make_prefix(start + timedelta(days=day), paths['NEMO'])
        path_NEMO_d = make_prefix(start + timedelta(days=day), paths['NEMO'],res='d')
        Ulist.append(path_NEMO + '_grid_U.nc')
        Vlist.append(path_NEMO + '_grid_V.nc')
        Wlist.append(path_NEMO + '_grid_W.nc')
        Tlist.append(path_NEMO + '_grid_T.nc')
        Waveslist.append(get_WW3_path(start + timedelta(days=day)))        

    # Load NEMO forcing 
    filenames = {
        'U': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Ulist},
        'V': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Vlist},
        'W': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Wlist},
        'Kz': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Wlist[0], 'data': Wlist},
        'T': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Tlist[0], 'data': Tlist},
        'S': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Tlist[0], 'data': Tlist},
        'ssh': {'lon': paths['coords'], 'lat': paths['coords'], 'data': Tlist},
        'R': {'lon': paths['coords'], 'lat': paths['coords'], 'depth': Tlist[0], 'data': Tlist},
        'Bathy' : {'lon': paths['coords'], 'lat': paths['coords'], 'data': paths['bat']},
        'gdepth' : {'lon': paths['coords'], 'lat': paths['coords'],'depth': Wlist[0], 'data': paths['mask']},
        'totdepth' : {'lon': paths['coords'], 'lat': paths['coords'], 'data': paths['mask']},
        'US' : {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist},
        'VS' : {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist},
        'WL' : {'lon': paths['coordsWW3'], 'lat': paths['coordsWW3'], 'data': Waveslist}
    }
    
    variables = {'U': 'vozocrtx', 'V': 'vomecrty','W': 'vovecrtz','T':'votemper','S':'vosaline','R':'sigma_theta',
        'US':'uuss','VS':'vuss','WL':'lm','Bathy':'Bathymetry','Kz':'vert_eddy_diff','ssh':'sossheig','totdepth':'totaldepth'}
        
    file_out,var_out = {},{}
    for var in varlist:
        file_out[var]=filenames[var]
        var_out[var]=variables[var]
    return file_out,var_out