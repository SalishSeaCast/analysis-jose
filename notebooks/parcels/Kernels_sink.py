def Buoyancy(particle, fieldset, time):
    '''Stokes law settling velocity'''
    if particle.sediment == 0:
        Rp = particle.ro-1000 #Density particle: LDPE (~920 kg/m3 ),PS (~150 kg/m3), PET (~1370 kg/m3).
        d = particle.diameter # particle diameter
        l = particle.length # particle length
        visc=1e-3 #average viscosity sea water
        z = particle.depth
        bath = 10*fieldset.mbathy[time, particle.depth, particle.lat, particle.lon]
        if  z > bath:
            particle.sediment = 1
        else:
            g = 9.8
            #ParcelsRandom.uniform(0,2)
            #t = fieldset.T[time, particle.depth, particle.lat, particle.lon]
            ro = fieldset.R[time, particle.depth, particle.lat, particle.lon]
            dro = Rp-ro
            #visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296)
            Ws= ((l/d)**-1.664)*0.079*((l**2)*g*(dro))/(visc)
            dz = Ws*particle.dt
        if dz+z > 0:
            particle.depth += dz 
        else:
            particle.depth = 0.5

def Stokes_drift(particle, fieldset, time):
    '''Stokes drift'''  
    if particle.sediment == 0:
        lat = particle.lat
        if lat > 48 and lat < 51: #Check that particle is inside WW3 data field
            Re = 6378137
            PI = math.pi
            z0 = particle.depth
            (us0, vs0, wl) = fieldset.stokes[time, particle.depth, particle.lat, particle.lon]
            k = (2*PI)/wl
            us = 180*(us0*exp(-abs(2*k*z0)))/(Re*PI*cos((particle.lat*PI)/180))
            vs = 180*(vs0*exp(-abs(2*k*z0)))/(Re*PI)
            particle.lon += us * particle.dt
            particle.lat += vs * particle.dt
        
def DeleteParticle(particle, fieldset, time):
    """Delete particle from OceanParcels simulation to avoid run failure
    """
    
    print(f'Particle {particle.id} lost !! [{particle.time}, {particle.depth}, {particle.lat}, {particle.lon}]')
    particle.delete()


def AdvectionRK4_3D(particle, fieldset, time):
    if particle.sediment == 0:
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
