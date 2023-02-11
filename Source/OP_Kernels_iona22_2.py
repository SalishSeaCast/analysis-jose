
def DeleteParticle(particle, fieldset, time):
    """Delete particle from OceanParcels simulation to avoid run failure"""
    print(f'Particle {particle.id} lost !! [{particle.time}, {particle.depth}, {particle.lat}, {particle.lon}]')
    particle.delete()
        
def Stokes_drift(particle, fieldset, time):
    """Stokes drift"""
    if particle.beached == 0:
        lat = particle.lat
        if lat > 48 and lat < 51: #Check that particle is inside WW3 data field
            deg2met_st = 111319.5
            latT_st = 0.6495 #cos(particle.lat*(math.pi/180))
            (us0, vs0, wl) = fieldset.stokes[time, particle.depth/dfactor, particle.lat, particle.lon]
            k = (2*math.pi)/wl
            us = (us0*exp(-math.fabs(2*k*particle.depth/dfactor)))/(deg2met_st*latT_st)
            vs = (vs0*exp(-math.fabs(2*k*particle.depth/dfactor)))/deg2met_st
            particle.lon += us * particle.dt 
            particle.lat += vs * particle.dt
        
def AdvectionRK4_3D(particle, fieldset, time):
    if particle.beached == 0: #Check particle is in the water column
        particle.tau += particle.dt
        if particle.tau > particle.dtmax:
            particle.delete()
        TD = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon]
        ssh = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon]
        dfactor = (1+ssh/TD)
        (u1, v1, w1) = fieldset.UVW[time, particle.depth/dfactor, particle.lat, particle.lon]
        lon1 = particle.lon + u1*.5*particle.dt
        lat1 = particle.lat + v1*.5*particle.dt
        dep1 = particle.depth + w1*.5*particle.dt
        (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1/dfactor, lat1, lon1]
        lon2 = particle.lon + u2*.5*particle.dt
        lat2 = particle.lat + v2*.5*particle.dt
        dep2 = particle.depth + w2*.5*particle.dt
        (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2/dfactor, lat2, lon2]
        lon3 = particle.lon + u3*particle.dt
        lat3 = particle.lat + v3*particle.dt
        dep3 = particle.depth + w3*particle.dt
        (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3/dfactor, lat3, lon3]
        particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.depth += (w1 + 2*w2 + 2*w3 + w4) / 6. * particle.dt


def turb_mix(particle,fieldset,time):
    """Vertical mixing and applying buoyancy"""
    if particle.beached==0:
        bath = fieldset.Bathymetria[time, particle.depth, particle.lat, particle.lon]
        if particle.depth/dfactor + 0.5 > bath: #Only calculate gradient of diffusion for particles deeper than 0.6 otherwise OP will check for particles outside the domain and remove it.
            Kzdz = 0
        else: 
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth/dfactor+0.5, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth/dfactor, particle.lat, particle.lon]) #forward difference 
        dgrad = Kzdz*particle.dt
        if particle.depth/dfactor+0.5*dgrad > 0.5:
            Kz = fieldset.vert_eddy_diff[time, particle.depth/dfactor+0.5*dgrad, particle.lat, particle.lon] #Vertical diffusivity SSC  #
        else:
            Kz = fieldset.vert_eddy_diff[time, 0.5, particle.lat, particle.lon] #Vertical diffusivity SSC  #

        Rr = ParcelsRandom.uniform(-1, 1)
        d_random = sqrt(3*2*Kz*particle.dt) * Rr
        dzs = dgrad + d_random
    ####Boundary conditions 
        Dlayer = 0.5*sqrt(Kz*particle.dt) #mixing layer dependant on Kz
        #First turbulent mixing second buoyancy.
        if dzs + particle.depth/dfactor > bath: #randomly in the water colum01n
            particle.depth = bath - Dlayer * ParcelsRandom.uniform(0, 1)
        elif particle.depth/dfactor + dzs < 0:
            particle.depth = Dlayer * ParcelsRandom.uniform(0, 1) #Well mixed boundary layer
        else:
            particle.depth += dzs #apply mixing
        if particle.depth/dfactor < 2:
            particle.beached = 6

def Beaching(particle, fieldset, time):
    """Horizontal mixing to impose beaching for particles reaching coast"""  
    if particle.beached == 0:
        kh = particle.Kh
        Rrx = ParcelsRandom.uniform(-1, 1)
        Rry = ParcelsRandom.uniform(-1, 1)
        d_x = sqrt(3*2*kh*particle.dt) * Rrx
        d_y = sqrt(3*2*kh*particle.dt) * Rry   
        deg2met = 111319.5
        latT = 0.6495
        d_randomx = particle.lon + d_x/(deg2met*latT)
        d_randomy = particle.lat + d_y/deg2met
        Sbh = fieldset.S[time, 1, d_randomy, d_randomx] #Check if particles reach coast (Salinity = 0)
        if particle.lat < 49.237 and particle.lon > -123.196 and particle.lat > 49.074:
            pass #Dont let particles beach inside the fraser river
        else:
            particle.lat = d_randomy
            particle.lon = d_randomx

def Unbeaching(particle, fieldset, time):
    '''Resuspension prob'''  
    if particle.beached == 1: 
        particle.tau += particle.dt
        if particle.tau > particle.dtmax:
            particle.delete()       
        Ub = particle.Ub*86400  #timescale unbeaching in seconds
        Pr = 1 - exp(-particle.dt/Ub)
        if ParcelsRandom.uniform(0,1)<Pr:
            particle.beached = 0
    elif particle.beached == 3:
        particle.tau += particle.dt