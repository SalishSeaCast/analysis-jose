def init_MP(particle, fieldset, time): 
    if particle.status==0: #Check if particle is initialized    
        '''Assign length and diameter to particles following input distribution''' 
        particle.status = 1 #free particle
        #Assign a value of diameter given size distribution
        d = ParcelsRandom.normalvariate(fieldset.diameter, fieldset.SDD)
        l = ParcelsRandom.normalvariate(fieldset.length, fieldset.SDL)
        #Do not allow negative dimensions
        while d <= fieldset.min_diameter:
            d = ParcelsRandom.normalvariate(fieldset.diameter, fieldset.SDD)
        while l <= fieldset.min_length:
            l = ParcelsRandom.normalvariate(fieldset.length, fieldset.SDL)
        particle.diameter = d
        particle.length = l

def Advection(particle, fieldset, time): 
    if particle.status == 1: #If particle is not active, do not move it
        ssh = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon,particle] #SSH(t)
        sshn = fieldset.sossheig[time+particle.dt, particle.depth, particle.lat, particle.lon,particle] #SSH(t+dt)
        td = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon,particle]#Total_depth
        particle.fact = (1+ssh/td) #VVl conversion factor
        VVL = (sshn-ssh)*particle.depth/td #VVl w velocity correction
        #4th order RK
        (u1, v1, w1) = fieldset.UVW[particle]
        lon1 = particle.lon + u1*.5*particle.dt
        lat1 = particle.lat + v1*.5*particle.dt
        dep1 = particle.depth + w1*.5*particle.dt/particle.fact
        (u2, v2, w2) = fieldset.UVW[time + .5 * particle.dt, dep1, lat1, lon1,particle]
        lon2 = particle.lon + u2*.5*particle.dt
        lat2 = particle.lat + v2*.5*particle.dt
        dep2 = particle.depth + w2*.5*particle.dt/particle.fact
        (u3, v3, w3) = fieldset.UVW[time + .5 * particle.dt, dep2, lat2, lon2,particle]
        lon3 = particle.lon + u3*particle.dt
        lat3 = particle.lat + v3*particle.dt
        dep3 = particle.depth + w3*particle.dt/particle.fact
        (u4, v4, w4) = fieldset.UVW[time + particle.dt, dep3, lat3, lon3,particle]
        particle_dlon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
        particle_dlat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
        particle_ddepth += (w1 + 2*w2 + 2*w3 + w4) /6. * particle.dt/particle.fact + VVL
        if particle_ddepth + particle.depth < 0:
            particle_ddepth = -(2*particle.depth + particle_ddepth) #reflection on surface

def vert_mix(particle,fieldset,time):
    """Vertical mixing"""
    if particle.status==1:
        #Vertical mixing
        if particle.depth + 0.5/particle.fact > td: #Only calculate gradient of diffusion for particles deeper than 0.5 otherwise OP will look for particles outside the domain raising an error.
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon,particle]-fieldset.vert_eddy_diff[time, particle.depth-0.5/particle.fact, particle.lat, particle.lon]) 
            #backwards difference 
        else: 
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth+0.5/particle.fact, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon,particle]) #forward difference 
        dgrad = Kzdz*particle.dt/particle.fact
        if particle.depth+(0.5*dgrad) > 0 and particle.depth+(0.5*dgrad) < td:
            Kz = fieldset.vert_eddy_diff[time, particle.depth+ 0.5*dgrad, particle.lat, particle.lon] 
            #Vertical diffusivity SSC  
        else:
            Kz = 0 #Otherwise it would check outside the domain

        Rr = ParcelsRandom.uniform(-1, 1)
        d_random = sqrt(3*2*Kz*particle.dt) * Rr/particle.fact
        dzs = (dgrad + d_random)
        
        #Apply turbulent mixing.
        if dzs + particle_ddepth + particle.depth > td: #crossed bottom in dt (Reflect)
            particle_ddepth = 2*td - (dzs + particle_ddepth + 2*particle.depth) #bounce on boundary/ no energy loss
        elif dzs + particle.depth + particle_ddepth < 0:#crossed surface in dt (#Well mixed boundary layer)
            particle_ddepth = -(dzs + 2*particle.depth+ particle_ddepth) 
        else:
            particle_ddepth += dzs #apply mixing  

def horz_mix(particle,fieldset,time):
    """Horizontal mixing (Beaching)"""
    if particle.status==1:
        Rrx = ParcelsRandom.uniform(-1, 1)
        Rry = ParcelsRandom.uniform(-1, 1)
        d_x = (sqrt(3*2*fieldset.kh*particle.dt) * Rrx)/fieldset.deglatT
        d_y = (sqrt(3*2*fieldset.kh*particle.dt) * Rry)/fieldset.deg2met   
        d_randomx = particle.lon + d_x + particle_dlon
        d_randomy = particle.lat + d_y + particle_dlat
        Sbh = 1
        Swh = 1
        if particle.depth < 5:
            Sbh = fieldset.sigma_theta[time, 0.51, d_randomy, d_randomx] #Check if particles reach coast at shallow water (SW Density = 0)
            if Sbh==0:
                particle.status = 2 #beached
        else:
            Swh = fieldset.sigma_theta[time, particle.depth, d_randomy, d_randomx]
            if Swh == 0:
                #Do not cross wall Keep particle in place
                particle.status = -1
                
        #Apply horizontal mixing (beaching for particles pushed through coast) 
        if particle.status == 1:
            particle_dlat += d_y
            particle_dlon += d_x
        elif particle.status == -1:
            particle.status = 1 #Do nothing, and resuspend. This particles crossed bathymetry at deep water

def Windage(particle, fieldset, time):
    if particle.status == 1 and particle.depth < 0.1:
        (uwind, vwind) = fieldset.wind[time, particle.depth, particle.lat, particle.lon+360,particle]
        particle_dlon += uwind * 0.03 * particle.dt 
        particle_dlat += vwind * 0.03 * particle.dt

def Stokes_drift(particle, fieldset, time):
    """Apply Stokes drift calculated by WW3"""
    if particle.status == 1 and particle.depth < 5:

        if particle.lat > 48 and particle.lat < 51: #particle is inside WW3 data field
            (us0, vs0, wl) = fieldset.stokes[time, 0, particle.lat, particle.lon,particle]
            k = (2*math.pi)/wl
            us = (us0*exp(-math.fabs(2*k*particle.depth)))
            vs = (vs0*exp(-math.fabs(2*k*particle.depth)))
        else:
            us = 0
            vs = 0
        particle_dlon += us * particle.dt
        particle_dlat += vs * particle.dt

def Buoyancy(particle, fieldset, time):
    """Calculating settling velocity using Komar cylinder or Fan super ellipsoids W_s"""
    if particle.status == 1: #Check particle is in the water column   
        deps = max(particle.depth,0.51) #Surface value of Tracer
        c = particle.dim1 # particle small semiaxis
        b = particle.dim2 # particle medium semiaxis
        a = particle.dim3 # particle large semiaxis
        e = particle.shape #particle roundness (0 is cube, 1 is sphere)
        rho_p = particle.density #particle density 

        g = 9.81 #gravity acceleration
        t = fieldset.votemper[time, deps, particle.lat, particle.lon] #Loading temperature from SSC
        rho_f = fieldset.sigma_theta[time, deps, particle.lat, particle.lon] + 1000 #Loading density sw from SSC
        mu = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296) #kinematic viscosity for Temp of SSC
        dr = math.fabs(rho_p-rho_f)

        #IS THE PARTICLE SINKING OR FLOATING
        sinking = 1
        if rho_p <= rho_f:
            sinking = -1 #raise ValueError

        EI = b / a 

        if EI <= 0.1: #CHECK if it's a fiber
                #fiber path through Komar
                l = a*2
                d = b + c
                W_s = ((l/d)**-1.664)*0.079*((l**2)*g*(dr))/(mu)
        else:     
            #fragment path through Fan
            FI = c / b
            # equivalent-volume sphere diameter
            Vol = 2.0*a*b*c * (e**2) * math.gamma(e/2 + 1)*math.gamma(e)/math.gamma(3*e/2 + 1) * math.gamma(e/2)*math.gamma(e/2)/math.gamma(e)
            d_eq = (6.0 * Vol / math.pi) ** (1.0 / 3.0)
    
            d_star = d_eq * (rho_f * dr * g / mu ** 2) ** (1.0 / 3.0)
    
            # sphere reference and shape correction
            u_star_S = 1.0 / (18.0 / d_star ** 2 + 2.412 / (4.0 * math.sqrt(d_star)))
            
            L = math.log10(d_star) #(Eq. 31b)
            C1 = 0.11 * L + 0.17
            C2 = 0.21 * L + 0.17
            C3 = 0.48 * L - 0.56
            C4 = 1.05 * L - 1.38
            
            u_star = EI ** C1 * FI ** C2 * (C3 * e ** 2 + C4 * e + 1.0 - C3 - C4) * u_star_S
            W_s = u_star * (dr * mu * g / rho_f ** 2) ** (1.0 / 3.0)
        
        dws = sinking * W_s*particle.dt/particle.fact
        
        #Apply buoyancy
        if particle.depth + dws + particle_ddepth > td:
            particle_ddepth = td-particle.depth #particle reaches the bottom 
            vtau_crit = (0.30 / (1.0 + 1.2 * d_star) + 0.055 * (1.0 - exp(-0.020 * d_star))) * g * d_eq
            vtau = fieldset.bottom_stress[time,particle.depth, particle.lat, particle.lon,particle]
            if vtau <= vtau_crit: #Assess critical shear here to avoid sedimentation
                particle.status = 3 #Trap particle in sediment

        elif particle.depth + dws + particle_ddepth < 0: #particle reaches the surface 
            particle_ddepth = - particle.depth 
        else:
            particle_ddepth += dws

def Unbeaching(particle, fieldset, time):
    '''Resuspension prob'''  
    if particle.status == 2:      
        if ParcelsRandom.uniform(0,1)<fieldset.PrUb:
            particle.status = 1

def resuspension(particle, fieldset, time):
    # proceed if particle is at the seabed
    if particle.status == 3:
        vtau_crit = fieldset.tau_crit #Value to check
        vtau = fieldset.Bot_stress[time,particle.depth, particle.lat, particle.lon,particle]
        if vtau <= vtau_crit:
            particle.status = 1
            #particle_ddepth += #!HOW MUCH should we resuspend?


def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        if Sbh == 0: #Ready to beach particles
            print('particle beached')
            particle.state = StatusCode.Success
            #particle.status = 2 #beached DO NOT MOVE outside domain
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

def export(particle,fieldset,time):
    if particle.status==1:#!need to add northern boundary too
        Pacific =  -particle.lat*0.84120957 -83.98027258 #Checking if particle gets too close to boundary JdF
        South = particle.lon*0.35157547 +90.26497859 #Checking southern boundary model
        if particle.lon<Pacific:
            particle.status = 5 #Exported to Pacific Ocean
        if particle.lat<South:
            if particle.depth < 5:
                particle.status =2 #Beached
            else:
                particle_dlat += 0.1*particle.dt/fieldset.deg2met 
                #Northward velocity to avoid particles from getting stuck on the wall only an issue at southernmost cell of SSC