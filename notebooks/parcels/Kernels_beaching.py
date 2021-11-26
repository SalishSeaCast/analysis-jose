def Buoyancy(particle, fieldset, time):
    '''Stokes law settling velocity'''
    Rp = particle.ro #Density particle: LDPE (~920 kg/m3 ),PS (~150 kg/m3), PET (~1370 kg/m3).
    ESD= particle.size #Size particle (ESD) equivalent spherical diameter
    visc=1e-3 #average viscosity sea water
     
    z = particle.depth
    u = fieldset.U[time, particle.depth, particle.lat, particle.lon]
    if  u == 0:
        dz = 0
    else:
        g = 9.8
        #t = fieldset.T[time, particle.depth, particle.lat, particle.lon]
        ro = fieldset.R[time, particle.depth, particle.lat, particle.lon]
        #visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296)
        Ws= ((ESD**2)*g*(Rp-ro))/(18*visc)
        WWS =(11.68 + 0.1991*ESD*1e6 + 0.0004*(ESD*1e6)**2- 0.0993*(Rp-ro) + 0.0002*(Rp-ro)**2)/86400
        if Ws > WWS:
            dz = WWS*particle.dt
        else:
            dz = Ws*particle.dt
    if dz+z > 0:
        particle.depth += dz 
    else:
        particle.depth = 0.5
        
def DeleteParticle(particle, fieldset, time):
    """Delete particle from OceanParcels simulation to avoid run failure
    """
    
    print(f'Particle {particle.id} lost !! [{particle.lon}, {particle.lat}, {particle.depth}, {particle.time}]')
    particle.delete()
    
def AdvectionRK4_3D(particle, fieldset, time):
    if particle.beached == 0:
        (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        lon1 = particle.lon + u1*.5*particle.dt
        lat1 = particle.lat + v1*.5*particle.dt
        dep1 = particle.depth + w1*.5*particle.dt
        (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
        lon2 = particle.lon + u2*.5*particle.dt
        lat2 = particle.lat + v2*.5*particle.dt
        dep2 = particle.depth + w2*.5*particle.dt
        (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
        lon3 = particle.lon + u3*particle.dt
        lat3 = particle.lat + v3*particle.dt
        dep3 = particle.depth + w3*particle.dt
        (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.depth += (w1 + 2*w2 + 2*w3 + w4) / 6. * particle.dt
        particle.beached = 2

def Stokes_drift(particle, fieldset, time):
    '''Stokes drift'''  
    lat = particle.lat
    if lat > 48 and lat < 51 and particle.beached == 0: #Check that particle is inside WW3 data field
        R = 6378137
        z = particle.depth
        us0 = fieldset.lm[time,particle.depth, particle.lat, particle.lon]
        vs0 = fieldset.vuss[time, particle.depth, particle.lat, particle.lon]
        wl = fieldset.lm[time, particle.depth, particle.lat, particle.lon]
        k = (2*math.pi)/wl
        us = 180*(us0*exp(-abs(2*k*z)))/(R*math.pi*cos((particle.lat*math.pi)/180))
        vs = 180*(vs0*exp(-abs(2*k*z)))/(R*math.pi)
        particle.lon += us * particle.dt
        particle.lat += vs * particle.dt
        particle.beached = 2

def BeachTesting_3D(particle, fieldset, time):
    if particle.beached == 2:
        (u, v, w) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        if fabs(u) < 1e-14 and fabs(v) < 1e-14:    
            particle.beached = 3
        else:
            particle.beached = 0