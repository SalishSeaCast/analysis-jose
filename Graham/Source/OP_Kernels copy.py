def Advection(particle, fieldset, time):
    if particle.status < 2: #Check particle isn't stuck
        if particle.status==0: #Check if particle is initialized    
            '''Assign length and diameter to particles following input distribution''' 
            particle.status = 1 #free particle
            SDD = 3.34e-6 #std Diameter
            SDL = 24.4e-5 #std length
            particle.diameter = ParcelsRandom.normalvariate(particle.diameter, SDD) #Randomly assign a value of diameter inside the Bamfield mesocosm size distribution
            particle.length = ParcelsRandom.normalvariate(particle.length, SDL) #Same for length
            '''Do not allow negative dimensions'''
            if particle.length < 0:
                particle.length = 4.5e-4-SDL
            if particle.diameter < 0:
                particle.length = 2.16e-5-SDD
        particle.tau += particle.dt #track age particle
        #if particle.tau > particle.dtmax: 
        #    '''If we set a maximum age for a particle it will be deleted here'''
        #    particle.delete()
        particle.dlon = 0 #initialize displacement lon
        particle.dlat = 0 #initialize displacement lat
        particle.ddepth = 0 #initialize displacement depth
        ssh = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon] #SSH(t)
        sshn = fieldset.sossheig[time+particle.dt, particle.depth, particle.lat, particle.lon] #SSH(t+dt)
        td = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon]#Total_depth
        particle.fact = (1+ssh/td)
        VVL = -((sshn-ssh)*particle.depth/(td))
        particle.cellvol = fieldset.volume[time, 0, particle.lat, particle.lon]
        (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
        lon1 = particle.lon + u1*.5*particle.dt
        lat1 = particle.lat + v1*.5*particle.dt
        dep1 = particle.depth + w1*.5*particle.dt/particle.fact
        (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
        lon2 = particle.lon + u2*.5*particle.dt
        lat2 = particle.lat + v2*.5*particle.dt
        dep2 = particle.depth + w2*.5*particle.dt/particle.fact
        (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
        lon3 = particle.lon + u3*particle.dt
        lat3 = particle.lat + v3*particle.dt
        dep3 = particle.depth + w3*particle.dt/particle.fact
        (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
        wa = (w1 + 2*w2 + 2*w3 + w4) /6.
        particle.wa = wa* particle.dt/particle.fact
        particle.dlon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle.dlat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle.ddepth += particle.wa - VVL

def Stokes_drift(particle, fieldset, time):
    """Apply Stokes drift calculated by WW3"""
    if particle.status == 1:
        lat = particle.lat
        if lat > 48 and lat < 51: #Check that particle is inside WW3 data field
            deg2met = 111319.5
            latT = 0.6495 #Average value for SoG #cos(particle.lat*(math.pi/180))
            (us0, vs0, wl) = fieldset.stokes[time, 0, particle.lat, particle.lon]
            k = (2*math.pi)/wl
            us = (us0*exp(-math.fabs(2*k*particle.depth)))/(deg2met*latT)
            vs = (vs0*exp(-math.fabs(2*k*particle.depth)))/deg2met
            particle.dlon += us * particle.dt
            particle.dlat += vs * particle.dt

def Buoyancy(particle, fieldset, time):
    """Calculating settling velocity using Komar cylinder Vs"""
    if particle.status == 1: #Check particle is in the water column   
        deps = max(particle.depth,0.51) #Surface value of Tracer node
        d = particle.diameter # particle diameter
        l = particle.length # particle length
        g = 9.81 #Gravity
        pro = 1350
        t = fieldset.votemper[time, deps, particle.lat, particle.lon] #Loading temperature from SSC
        ro = fieldset.sigma_theta[time, deps, particle.lat, particle.lon] #Loading density sw from SSC
        visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296) #kinematic viscosity for Temp of SSC
        Ws= ((l/d)**-1.664)*0.079*((l**2)*g*(pro-1000-ro))/(visc)
        #print('settling velocity in m/h =',str(Ws*3600))
        particle.ws = Ws
        dws = Ws*particle.dt/particle.fact

def turb_mix(particle,fieldset,time):
    """Vertical mixing"""
    if particle.status==1:
        #Vertical mixing
        if particle.depth + 0.5/particle.fact > td: #Only calculate gradient of diffusion for particles deeper than 0.5 otherwise OP will check for particles outside the domain and remove it.
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth-0.5/particle.fact, particle.lat, particle.lon]) #backwards difference 
        else: 
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth+0.5/particle.fact, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]) #forward difference 
        dgrad = Kzdz*particle.dt
        if particle.depth+0.5*dgrad > 0 and particle.depth+0.5*dgrad < td:
            Kz = fieldset.vert_eddy_diff[time, particle.depth+ 0.5*dgrad/particle.fact, particle.lat, particle.lon] #Vertical diffusivity SSC  
        else:
            Kz = fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon] 
            #print('out of bounds mixing=',Kz)

        Rr = ParcelsRandom.uniform(-1, 1)
        d_random = sqrt(3*2*Kz*particle.dt) * Rr
        dzs = (dgrad + d_random)/particle.fact
        #print('mixing displacement ==', str(dzs*3600))
        #Horizontal mixing (Beaching BC)
        kh = 1.5 #SSC constant value
        Rrx = ParcelsRandom.uniform(-1, 1)
        Rry = ParcelsRandom.uniform(-1, 1)
        d_x = (sqrt(3*2*kh*particle.dt) * Rrx)/(deg2met*latT)
        d_y = (sqrt(3*2*kh*particle.dt) * Rry)/deg2met   
        d_randomx = particle.lon + d_x + particle.dlon
        d_randomy = particle.lat + d_y + particle.dlat
        Sbh = fieldset.vosaline[time, 0.51, d_randomy, d_randomx] #Check if particles reach coast at surface (Salinity = 0)
    
def Displacement(particle,fieldset,time):
    ''''Apply movement calculated by other kernels'''
    if particle.status == 1:
        #Apply buoyancy
        if particle.depth + dws + particle.ddepth > td:
            particle.ddepth = td-particle.depth #particle reaches the bottom 
            particle.status = 3 #Trap particle in sediment (sticky bottom)
            if ParcelsRandom.uniform(0,1) > particle.alpha:
                particle.status = 1
        elif particle.depth + dws + particle.ddepth < 0: #particle reaches the surface 
            particle.ddepth = -particle.depth #! check if loosing too many particles through surface
        else:
            particle.ddepth += dws

        #Apply turbulent mixing.
        if dzs + particle.ddepth + particle.depth > td: #crossed bottom in dt (Reflect)
            particle.ddepth = 2*td - (dzs + particle.ddepth + 2*particle.depth) #bounce on boundary/ no energy loss
        elif dzs + particle.depth < 0:
            particle.ddepth = -(dzs + 2*particle.depth) #Well mixed boundary layer
        else:
            particle.ddepth += dzs #apply mixing  
        #Apply horizontal mixing (beaching for particles pushed through coast) 
        if Sbh == 0:
            particle.status = 2 #beached DO NOT MOVE outside domain
        else:
            particle.dlat += d_y
            particle.dlon += d_x
    
        #UPDATE LOCATION   
        particle.depth += particle.ddepth
        particle.lat += particle.dlat
        particle.lon += particle.dlon
        particle.time += particle.dt


def Unbeaching(particle, fieldset, time):
    '''Resuspension prob'''  
    if particle.status == 2: 
        particle.tau += particle.dt
        # if particle.tau > particle.dtmax:
        #     particle.delete()       
        Ub = particle.Ub*86400  #timescale unbeaching in seconds
        Pr = 1 - exp(-particle.dt/Ub)
        if ParcelsRandom.uniform(0,1)<Pr:
            particle.status = 1
    elif particle.status > 2:
        particle.tau += particle.dt


def DeleteParticle(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
       particle.delete()