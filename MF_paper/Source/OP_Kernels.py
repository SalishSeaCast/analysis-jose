def Advection(particle, fieldset, time): 
    if particle.status < 2: #Check particle isn't stuck
        if particle.status==0: #Check if particle is initialized    
            '''Assign length and diameter to particles following input distribution''' 
            particle.status = 1 #free particle
            SDD = 3.34e-6 #std Diameter
            SDL = 24.4e-5 #std length
            #Assign a value of diameter given size distribution
            particle.diameter = ParcelsRandom.normalvariate(particle.diameter, SDD) 
            particle.length = ParcelsRandom.normalvariate(particle.length, SDL)
            '''Do not allow negative dimensions'''
            if particle.length < 0:
                particle.length = 4.5e-4-SDL
            if particle.diameter < 0:
                particle.length = 2.16e-5-SDD
        #SOME CONSTANTS NEEDED LATER 
        deg2met = 111319.5
        latT = 0.6495 #Average value for SoG #cos(particle.lat*(math.pi/180))
        ssh = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon] #SSH(t)
        sshn = fieldset.sossheig[time+particle.dt, particle.depth, particle.lat, particle.lon] #SSH(t+dt)
        td = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon]#Total_depth
        particle.fact = (1+ssh/td) #VVl conversion factor
        VVL = (sshn-ssh)*particle.depth/td #VVl w velocity correction
        #4th order RK
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
        particle_dlon = (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle_dlat = (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle_ddepth = particle.wa + VVL
        if particle_ddepth + particle.depth < 0:
            particle_ddepth = -(2*particle.depth + particle_ddepth) #reflection on surface

def Stokes_drift(particle, fieldset, time):
    """Apply Stokes drift calculated by WW3"""
    if particle.status == 1:
        if particle.depth < 5:
            if particle.lat > 48 and particle.lat < 51: #particle is inside WW3 data field
                (us0, vs0, wl) = fieldset.stokes[time, 0, particle.lat, particle.lon]
                k = (2*math.pi)/wl
                us = (us0*exp(-math.fabs(2*k*particle.depth)))/(deg2met*latT)
                vs = (vs0*exp(-math.fabs(2*k*particle.depth)))/deg2met
            else:
                us = 0
                vs = 0
        else:
            us = 0
            vs = 0
        particle_dlon += us * particle.dt
        particle_dlat += vs * particle.dt
    

def Buoyancy(particle, fieldset, time):
    """Calculating settling velocity using Komar cylinder Ws"""
    if particle.status == 1: #Check particle is in the water column   
        deps = max(particle.depth,0.51) #Surface value of Tracer node
        d = particle.diameter # particle diameter
        l = particle.length # particle length
        g = 9.81 #Gravity
        pro = 1350 #average value density PET
        t = fieldset.votemper[time, deps, particle.lat, particle.lon] #Loading temperature from SSC
        ro = fieldset.sigma_theta[time, deps, particle.lat, particle.lon] #Loading density sw from SSC
        visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296) #kinematic viscosity for Temp of SSC
        Ws= ((l/d)**-1.664)*0.079*((l**2)*g*(pro-1000-ro))/(visc)
        dws = Ws*particle.dt/particle.fact
        particle.ws = dws

def export(particle,fieldset,time):
    if particle.status==1:
        test =  -particle.lat*0.84120957 -83.98027258 #Checking if particle gets too close to boundary JdF
        test2 = particle.lon*0.35157547 +90.26497859 #Checking southern boundary model
        if particle.lon<test:
            particle.status = 5 #Exported to Pacific Ocean
        if particle.lat<test2:
            if particle.depth < 5:
                particle.status =2 #Beached
            else:
                particle_dlat += 0.1*particle.dt/deg2met 
                #Northward velocity to avoid particles from getting stuck on the wall only an issue at southernmost cell of SSC
        

def turb_mix(particle,fieldset,time):
    """Vertical mixing"""
    if particle.status==1:
        #Vertical mixing
        if particle.depth + 0.5/particle.fact > td: #Only calculate gradient of diffusion for particles deeper than 0.5 otherwise OP will look for particles outside the domain raising an error.
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth-0.5/particle.fact, particle.lat, particle.lon]) 
            #backwards difference 
        else: 
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth+0.5/particle.fact, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]) #forward difference 
        dgrad = Kzdz*particle.dt/particle.fact
        if particle.depth+(0.5*dgrad) > 0 and particle.depth+(0.5*dgrad) < td:
            Kz = fieldset.vert_eddy_diff[time, particle.depth+ 0.5*dgrad, particle.lat, particle.lon] 
            #Vertical diffusivity SSC  
        else:
            Kz = 0 #Otherwise it would check outside the domain

        Rr = ParcelsRandom.uniform(-1, 1)
        d_random = sqrt(3*2*Kz*particle.dt) * Rr/particle.fact
        dzs = (dgrad + d_random)
        particle.wm = dzs

        #Horizontal mixing (Beaching BC)
        kh = 1.5 #SSC constant value
        Rrx = ParcelsRandom.uniform(-1, 1)
        Rry = ParcelsRandom.uniform(-1, 1)
        d_x = (sqrt(3*2*kh*particle.dt) * Rrx)/(deg2met*latT)
        d_y = (sqrt(3*2*kh*particle.dt) * Rry)/deg2met   
        d_randomx = particle.lon + d_x + particle_dlon
        d_randomy = particle.lat + d_y + particle_dlat
        Sbh = 1
        Swh = 1
        if particle.depth < 5:
            Sbh = fieldset.sigma_theta[time, 0.51, d_randomy, d_randomx] #Check if particles reach coast at shallow water (SW Density = 0)
            if Sbh==0:
                particle.status = 2 
        else:
            Swh = fieldset.sigma_theta[time, particle.depth, d_randomy, d_randomx]
            if Swh == 0:
                #Do not cross wall Keep particle in place
                particle.status = -1

    
def Displacement(particle,fieldset,time):
    ''''Apply movement calculated by other kernels'''
    if particle.status == 1:
        #Apply buoyancy
        if particle.depth + dws + particle_ddepth > td:
            particle_ddepth = td-particle.depth #particle reaches the bottom 
            particle.status = 3 #Trap particle in sediment (sticky bottom)
            if ParcelsRandom.uniform(0,1) > particle.alpha:
                particle.status = 1 #resuspension don't allow particle to stick
        elif particle.depth + dws + particle_ddepth < 0: #particle reaches the surface 
            particle_ddepth = - particle.depth 
        else:
            particle_ddepth += dws

        #Apply turbulent mixing.
        if dzs + particle_ddepth + particle.depth > td: #crossed bottom in dt (Reflect)
            particle_ddepth = 2*td - (dzs + particle_ddepth + 2*particle.depth) #bounce on boundary/ no energy loss
        elif dzs + particle.depth + particle_ddepth < 0:
            particle_ddepth = -(dzs + 2*particle.depth+ particle_ddepth) #Well mixed boundary layer
        else:
            particle_ddepth += dzs #apply mixing  
        #Apply horizontal mixing (beaching for particles pushed through coast) 
        if particle.status == 1:
            particle_dlat += d_y
            particle_dlon += d_x
        elif particle.status == -1:
            particle.status = 1 #Do nothing, and resuspend. This particles crossed bathymetry at deep water
        

def Unbeaching(particle, fieldset, time):
    '''Resuspension prob'''  
    if particle.status == 2:      
        Ub = particle.Ub * 86400  #timescale unbeaching in seconds
        Pr = 1 - exp(-particle.dt/Ub)
        if ParcelsRandom.uniform(0,1)<Pr:
            particle.status = 1


def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        if Sbh == 0: #Ready to beach particles
            print('particle beached')
            particle.state = StatusCode.Success
            particle.status = 2 #beached DO NOT MOVE outside domain
        elif Swh ==0:
            print('particle hit wall') #Do not cross wall Keep particle in place
            particle.state = StatusCode.Success
        else:
            print('lost Particle')
            particle.status = 10 #lost
            particle.delete()
        
def KeepInOcean(particle, fieldset, time):
    if particle.state == StatusCode.ErrorThroughSurface:
        particle.depth = 0.0
        particle.state = StatusCode.Success
