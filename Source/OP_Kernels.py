def Advection(particle, fieldset, time):
    if particle.beached == 0: #Check particle is in the water column
        if particle.tau==0: #Check age particle is 0    
            '''Assign length and diameter to particles following input distribution''' 
            particle.diameter = ParcelsRandom.normalvariate(particle.diameter, particle.SDD) #Randomly assign a value of diameter inside the Bamfield mesocosm size dist
            particle.length = ParcelsRandom.normalvariate(particle.length, particle.SDL) #Same for length
            '''Check for negative dimensions'''
            if particle.length < 0:
                particle.length = 4.5e-4-particle.SDL
            if particle.diameter < 0:
                particle.length = 2.16e-5-particle.SDD
        particle.tau += particle.dt
        if particle.tau > particle.dtmax: #If not specified default value set to 14 years in OP_functions
            '''If we set a maximum age for a particle it will be deleted here'''
            particle.delete()
        if particle.surf == 1: #If particle is floating no effect of w (reduce computation) RK4_2D
            particle.depth = 0
            (u1, v1) = fieldset.UV[time, particle.depth, particle.lat, particle.lon]
            lon1 = particle.lon + u1*.5*particle.dt
            lat1 = particle.lat + v1*.5*particle.dt
            (u2, v2) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat1, lon1]
            lon2 = particle.lon + u2*.5*particle.dt
            lat2 = particle.lat + v2*.5*particle.dt
            (u3, v3) = fieldset.UV[time + .5 * particle.dt, particle.depth, lat2, lon2]
            lon3 = particle.lon + u3*particle.dt
            lat3 = particle.lat + v3*particle.dt
            (u4, v4) = fieldset.UV[time + particle.dt, particle.depth, lat3, lon3]
            particle.wa = 0
            particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
            particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        else: #RK4 3d + VVL correction
            ssh = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon] #SSH(t)
            sshn = fieldset.sossheig[time+particle.dt, particle.depth, particle.lat, particle.lon] #SSH(t+dt)
            td = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon]#Total_depth
            factor = (1+ssh/td)
            VVL = ((sshn-ssh)*particle.depth/(td+ssh))
            particle.fact = factor
            particle.cellvol = fieldset.volume[time, 0, particle.lat, particle.lon]
            (u1, v1, w1) = fieldset.UVW[time, particle.depth, particle.lat, particle.lon]
            lon1 = particle.lon + u1*.5*particle.dt
            lat1 = particle.lat + v1*.5*particle.dt
            dep1 = particle.depth + w1*.5*particle.dt/factor
            (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1]
            lon2 = particle.lon + u2*.5*particle.dt
            lat2 = particle.lat + v2*.5*particle.dt
            dep2 = particle.depth + w2*.5*particle.dt/factor
            (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2]
            lon3 = particle.lon + u3*particle.dt
            lat3 = particle.lat + v3*particle.dt
            dep3 = particle.depth + w3*particle.dt/factor
            (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3]
            wa = (w1 + 2*w2 + 2*w3 + w4) /6.
            particle.wa = wa
            particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
            particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
            particle.depth += wa * particle.dt/factor - VVL

def Stokes_drift(particle, fieldset, time):
    """Apply Stokes drift calculated by WW3"""
    if particle.beached == 0:
        lat = particle.lat
        if lat > 48 and lat < 51: #Check that particle is inside WW3 data field
            deg2met = 111319.5
            latT = 0.6495 #Average value for SoG #cos(particle.lat*(math.pi/180))
            (us0, vs0, wl) = fieldset.stokes[time, 0, particle.lat, particle.lon]
            k = (2*math.pi)/wl
            us = (us0*exp(-math.fabs(2*k*particle.depth)))/(deg2met*latT)
            vs = (vs0*exp(-math.fabs(2*k*particle.depth)))/deg2met
            particle.lon += us * particle.dt 
            particle.lat += vs * particle.dt

def Buoyancy(particle, fieldset, time):
    """Calculating settling velocity using Komar cylinder Vs"""
    if particle.beached == 0: #Check particle is in the water column   
        deps = max(particle.depth,0.51)
        bath = fieldset.Bathymetry[time, 0, particle.lat, particle.lon]  
        d = particle.diameter # particle diameter
        l = particle.length # particle length
        g = 9.8 #Gravity
        t = fieldset.votemper[time, deps, particle.lat, particle.lon] #Loading temperature from SSC
        ro = fieldset.sigma_theta[time, deps, particle.lat, particle.lon] #Loading density sw from SSC
        visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296) #kinematic viscosity for Temp of SSC
        Ws= ((l/d)**-1.664)*0.079*((l**2)*g*(particle.ro-1000-ro))/(visc)
        particle.ws = Ws/factor    

def turb_mix(particle,fieldset,time):
    """Vertical mixing"""
    if particle.beached==0:
        #Vertical mixing
        if particle.depth + 0.5/factor > bath: #Only calculate gradient of diffusion for particles deeper than 0.6 otherwise OP will check for particles outside the domain and remove it.
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth-0.5/factor, particle.lat, particle.lon]) #backwards difference 
        else: 
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth+0.5/factor, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]) #forward difference 
        dgrad = Kzdz*particle.dt
        if particle.depth+0.5*dgrad > 0 and particle.depth+0.5*dgrad < bath:
            Kz = fieldset.vert_eddy_diff[time, particle.depth+ 0.5*dgrad/factor, particle.lat, particle.lon] #Vertical diffusivity SSC  
        else:
            Kz = fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon] 

        Rr = ParcelsRandom.uniform(-1, 1)
        d_random = sqrt(3*2*Kz*particle.dt) * Rr
        dzs = dgrad + d_random
        Dlayer = 0.5*sqrt(Kz*particle.dt) #mixing layer dependant on Kz
        #Horizontal mixing (Beaching BC)
        kh = particle.Kh
        Rrx = ParcelsRandom.uniform(-1, 1)
        Rry = ParcelsRandom.uniform(-1, 1)
        d_x = sqrt(3*2*kh*particle.dt) * Rrx
        d_y = sqrt(3*2*kh*particle.dt) * Rry   
        d_randomx = particle.lon + d_x/(deg2met*latT)
        d_randomy = particle.lat + d_y/deg2met
        Sbh = fieldset.vosaline[time, 1, d_randomy, d_randomx] #Check if particles reach coast at surface (Salinity = 0)
    
def Displacement(particle,fieldset,time):
    ''''Apply movement calculated by other kernels'''
    if particle.beached==0 and particle.surf == 0:
        #Apply turbulent mixing.
        if dzs/factor + particle.depth > bath: #randomly in the water column
            particle.depth = bath - (Dlayer/factor) * ParcelsRandom.uniform(0, 1) #Well mixed boundary layer
        elif dzs/factor + particle.depth < 0:
            particle.depth = (Dlayer/factor) * ParcelsRandom.uniform(0, 1) #Well mixed boundary layer
        else:
            particle.depth += dzs/factor #apply mixing  
        #Apply horizontal mixing (beaching for particles pushed through coast) 
        if particle.lat < 49.237 and particle.lon > -123.196 and particle.lat > 49.074:
            pass #Dont let particles beach inside the Fraser river
        elif Sbh == 0:
            particle.beached = 1
        else:
            particle.lat = d_randomy
            particle.lon = d_randomx
        #Apply buoyancy
        if particle.depth + Ws*particle.dt/factor > bath:
            particle.beached = 3 #Trap particle in sediment (sticky bottom)
        elif particle.depth + Ws*particle.dt/factor < 0: 
            particle.depth = 0 #! check if loosing too many particles through surface
        else:
            particle.depth += Ws*particle.dt/factor
        


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
    elif particle.beached == 0:
        pass
    else:
        particle.tau += particle.dt

def DeleteParticle(particle, fieldset, time):
    """Delete particle from OceanParcels simulation to avoid run failure"""
    print(f'Particle {particle.id} lost !! [{particle.time}, {particle.depth}, {particle.lat}, {particle.lon}]')
    particle.beached = 6
