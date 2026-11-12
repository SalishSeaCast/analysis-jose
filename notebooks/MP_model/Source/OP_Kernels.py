def vert_mix(particle,fieldset,time):
    """Vertical mixing runs every 5 seconds to solve kzdz"""
    if particle.status==1:
        ssh = fieldset.sossheig[time, particle.depth, particle.lat, particle.lon,particle] #SSH(t)
        td = fieldset.totaldepth[time, particle.depth, particle.lat, particle.lon,particle]#Total_depth
        particle.fact = (1+ssh/td) #VVl conversion factor

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
            Kz = 0.0 #Otherwise it would check outside the domain

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

def Advection(particle, fieldset, time): 
    if particle.status == 1 and time - particle.last_update_dyn >= fieldset.dyn_dt: #If particle is not active, do not move it
        #4th order RK
        (u1, v1, w1) = fieldset.UVW[particle]
        lon1 = particle.lon + u1*.5*fieldset.dyn_dt
        lat1 = particle.lat + v1*.5*fieldset.dyn_dt
        dep1 = particle.depth + w1*.5*fieldset.dyn_dt/particle.fact
        (u2, v2, w2) = fieldset.UVW[time + .5 * fieldset.dyn_dt, dep1, lat1, lon1,particle]
        lon2 = particle.lon + u2*.5*fieldset.dyn_dt
        lat2 = particle.lat + v2*.5*fieldset.dyn_dt
        dep2 = particle.depth + w2*.5*fieldset.dyn_dt/particle.fact
        (u3, v3, w3) = fieldset.UVW[time + .5 * fieldset.dyn_dt, dep2, lat2, lon2,particle]
        lon3 = particle.lon + u3*fieldset.dyn_dt
        lat3 = particle.lat + v3*fieldset.dyn_dt
        dep3 = particle.depth + w3*fieldset.dyn_dt/particle.fact
        (u4, v4, w4) = fieldset.UVW[time + fieldset.dyn_dt, dep3, lat3, lon3,particle]
        particle_dlon += (u1 + 2*u2 + 2*u3 + u4) / 6. * fieldset.dyn_dt
        particle_dlat += (v1 + 2*v2 + 2*v3 + v4) / 6. * fieldset.dyn_dt
        particle_ddepth += (w1 + 2*w2 + 2*w3 + w4) / 6. * fieldset.dyn_dt/particle.fact 
        if particle_ddepth + particle.depth < 0:
            particle_ddepth = -(2*particle.depth + particle_ddepth) #reflection on surface

def horz_mix(particle,fieldset,time):
    """Horizontal mixing (Beaching)"""
    if particle.status==1 and time - particle.last_update_dyn >= fieldset.dyn_dt:
        Rrx = ParcelsRandom.uniform(-1, 1)
        Rry = ParcelsRandom.uniform(-1, 1)
        particle_dlon += (math.sqrt(3*2*fieldset.kh*fieldset.dyn_dt) * Rrx) / fieldset.deglatT
        particle_dlat += (math.sqrt(3*2*fieldset.kh*fieldset.dyn_dt) * Rry) / fieldset.deg2met


def LandCheck(particle, fieldset, time):
    """Cancel any horizontal move whose destination is a land cell (td == 0)."""
    if particle.status == 1 and time - particle.last_update_dyn >= fieldset.dyn_dt:
        td2 = fieldset.totaldepth[time, particle.depth, particle.lat + particle_dlat, particle.lon + particle_dlon, particle]
        if td2 == 0.0:
            particle_dlon = 0.0          # stay in the last wet position
            particle_dlat = 0.0
            if particle.depth < 5:
                particle.status = 2      # beached
        
def Stokes_Windage(particle, fieldset, time):
    """Apply Stokes drift calculated by WW3 and HRDPS Wind drag"""
    if particle.status == 1 and particle.depth < 5 and time - particle.last_update_dyn >= fieldset.dyn_dt:
        us = 0.0
        vs = 0.0
        if particle.lat > 48. and particle.lat < 51.: #particle is inside WW3 data field
            (us0, vs0, wl) = fieldset.stokes[time, 0, particle.lat, particle.lon,particle]
            if wl>0:
                k = (2*math.pi)/wl
                decay = math.exp(-math.fabs(2*k*particle.depth))
                us = us0*decay
                vs = vs0*decay
        particle_dlon += us * fieldset.dyn_dt
        particle_dlat += vs * fieldset.dyn_dt

        if particle.depth < 0.1: #Windage
            (uwind, vwind) = fieldset.wind[time, 0, particle.lat, particle.lon+360,particle]
            particle_dlon += uwind * 0.03 * fieldset.dyn_dt
            particle_dlat += vwind * 0.03 * fieldset.dyn_dt
            

def Biofouling_ws(particle, fieldset, time):
    if particle.status == 1 and time - particle.last_update_biol >= fieldset.biol_dt:
        HBAC = fieldset.heterotrophic_bacteria[time, particle.depth, particle.lat, particle.lon, particle]
        PHYT = fieldset.diatoms[time, particle.depth, particle.lat, particle.lon, particle]
        BP = fieldset.HetHBAC[time, particle.depth, particle.lat, particle.lon, particle]
        PP = fieldset.PPDIAT[time, particle.depth, particle.lat, particle.lon, particle]
        M_hbac = fieldset.MORTHBAC[time, particle.depth, particle.lat, particle.lon, particle]
        M_ph = fieldset.MORTDIAT[time, particle.depth, particle.lat, particle.lon, particle]

        growth_hb = 0.0
        growth_ph = 0.0
        mortality_hb = 0.0
        mortality_ph = 0.0

        if HBAC > 0.0:
            growth_hb = ((BP/HBAC)/86400.) * particle.N_hb
            mortality_hb = ((M_hbac/HBAC)/86400.) * particle.N_hb

        if PHYT > 0.0:
            growth_ph = ((PP/PHYT)/86400.) * particle.N_ph
            mortality_ph = ((M_ph/PHYT)/86400.) * particle.N_ph

        C_hb = HBAC * fieldset.N_hb_cell 
        C_ph = PHYT * fieldset.N_ph_cell 

        c = particle.dim1 # particle small semiaxis
        b = particle.dim2 # particle medium semiaxis
        a = particle.dim3 # particle large semiaxis
        e = particle.shape
        fe = e**2 * math.exp(3*math.lgamma(e/2) - math.lgamma(3*e/2)) #gamma == lgamma

        r_eq = (((4.*(a+particle.bioth)*(b+particle.bioth)*(c+particle.bioth)*fe)/math.pi) ** (1./3.))/2.
        r_hb = fieldset.r_hbac #!find representative values
        r_ph = fieldset.r_phyt #!find representative values

        deps = max(particle.depth,0.51) #Surface value of Tracer
        t = fieldset.votemper[time, deps, particle.lat, particle.lon] #temperature from SSC
        mu = 4.2844e-5 + 1./(0.157*((t + 64.993)**2)-91.296) #kinematic muosity for Temp of SSC
        gamma_s = 2. # Shear rate [s-1] #!sensitivity
        coll_eff = 1. #collision efficiency #!sensitivity

        # Relative settling velocities [m s-1]
        vs_hb = math.fabs(particle.ws - fieldset.vs_hbac)#!find representative values
        vs_ph = math.fabs(particle.ws - fieldset.vs_diat)#!find representative values

        # collision rate
        D_pl = fieldset.kB * (t + 273.15) / (6. * math.pi * mu * r_eq)
        D_hb = fieldset.kB * (t + 273.15) / (6. * math.pi * mu * r_hb)
        D_ph = fieldset.kB * (t + 273.15) / (6. * math.pi * mu * r_ph)

        beta_brown_hb = 4. * math.pi * (D_pl + D_hb) * (r_eq + r_hb)
        beta_brown_ph = 4. * math.pi * (D_pl + D_ph) * (r_eq + r_ph)
        beta_shear_hb = 1.3 * gamma_s * (r_eq + r_hb)**3
        beta_shear_ph = 1.3 * gamma_s * (r_eq + r_ph)**3
        beta_set_hb = 0.5 * math.pi * r_eq**2 * vs_hb
        beta_set_ph = 0.5 * math.pi * r_eq**2 * vs_ph

        beta_hb = beta_brown_hb + beta_shear_hb + beta_set_hb
        beta_ph = beta_brown_ph + beta_shear_ph + beta_set_ph

        R_coll_hb = beta_hb * C_hb * coll_eff 
        R_coll_ph = beta_ph * C_ph * coll_eff 

        # Detachment [cells s-1]
        detachment_hb = 0. #!sensitivity
        detachment_ph = 0.

        # Update attached populations
        particle.N_hb += (R_coll_hb + growth_hb - mortality_hb - detachment_hb)*fieldset.biol_dt
        particle.N_ph += (R_coll_ph + growth_ph - mortality_ph - detachment_ph)*fieldset.biol_dt

        # Prevent negative populations
        if particle.N_hb < 0.:
            particle.N_hb = 0.

        if particle.N_ph < 0.:
            particle.N_ph = 0.

        # fouling thickness (Newton on t^3 + A t^2 + B t = K)
        Vol_hbac = fieldset.V_hb * particle.N_hb
        Vol_phyto = fieldset.V_ph * particle.N_ph
        V_bio = Vol_hbac + Vol_phyto
        A = a + b + c
        B = a*c + b*c + a*b
        K = 3.*V_bio/(2.*fe)
        th = 0.
        if K > 0.:
            th = K/B
            i = 0
            while i < 6:
                th -= (th*(th*(th + A) + B) - K)/(3*th*th + 2*A*th + B)
                i += 1
        particle.bioth = th

        #Now let's calculate settling velocity with updated dimensions
        cf = c + th
        bf = b + th
        af = a + th
        rho_p = (particle.density*particle.Vol + fieldset.rho_hb*Vol_hbac + fieldset.rho_ph*Vol_phyto)/ (particle.Vol+ V_bio)
        g = fieldset.gravity 
        rho_f = fieldset.sigma_theta[time, deps, particle.lat, particle.lon] + 1000 #Loading density sw from SSC
        dr = math.fabs(rho_p-rho_f)

        if dr < 1e-6:
            particle.ws = 0.0
            particle.vtau_crit = 0.0
        else:
        #IS THE PARTICLE SINKING OR FLOATING
            sinking = 1
            if rho_p <= rho_f:
                sinking = -1 
        
            EI = bf/af 
            FI = cf/bf
            d_eq = ((4*af*bf*cf*fe)/math.pi)**(1./3.)
            d_star = d_eq * (rho_f*dr*g /mu**2)**(1./3.)
            L = math.log10(d_star) #(Eq. 31b)
            C1 = 0.11 * L + 0.17
            C2 = 0.21 * L + 0.17
            C3 = -0.48 * L + 0.56
            C4 = 1.05 * L - 0.38
            
            u_star_S = 1/(18/d_star**2 + (2.412/(4*math.sqrt(d_star))))
            u_star = EI**C1 * FI**C2 * (C3*e**2 + C4*e + 1 - C3 - C4) * u_star_S
            particle.ws = sinking * u_star * (dr*mu*g / rho_f**2)**(1./3.)
            particle.vtau_crit = (0.30 / (1.0 + 1.2 * d_star) + 0.055 * (1.0 - exp(-0.020 * d_star))) * g * d_eq * dr 
        particle.last_update_biol = time

def Buoyancy(particle, fieldset, time):
    """Calculating settling velocity using super ellipsoids W_s"""
    if particle.status == 1 and time - particle.last_update_dyn >= fieldset.dyn_dt: #Only update w_s once every hour: #Check particle is in the water column       
        dws = particle.ws*fieldset.dyn_dt/particle.fact
        #Apply buoyancy
        if particle.depth + dws + particle_ddepth > td:
            particle_ddepth = td-particle.depth #particle reaches the bottom 
            if particle.ws > 0:   # only sinking particles can settle into the sediment
                vtau = fieldset.bottom_stress[time,particle.depth, particle.lat, particle.lon,particle]
                if vtau <= particle.vtau_crit: #Assess critical shear here to avoid sedimentation
                    particle.status = 3 #Trap particle in sediment

        elif particle.depth + dws + particle_ddepth < 0: #particle reaches the surface 
            particle_ddepth = - particle.depth 
        else:
            particle_ddepth += dws
        particle.last_update_dyn = time

def Unbeaching(particle, fieldset, time):
    '''Resuspension prob'''  
    if particle.status == 2: #todo if want to increase dt we need to calculate PrUb accordingly
        if ParcelsRandom.uniform(0., 1.)<fieldset.PrUb:
            particle.status = 1

def resuspension(particle, fieldset, time):
    if particle.status == 3 and time - particle.last_update_res >= fieldset.dyn_dt:
        vtau = fieldset.bottom_stress[time, particle.depth, particle.lat, particle.lon, particle]
        if vtau > particle.vtau_crit:  #todo erode_factor set to 1.0 consider increasing to reduce resuspension
            particle.status = 1
            lift = -fieldset.resusp_h * math.log(1. - ParcelsRandom.uniform(0., 1.))   #!sensitivity 0.1, 0.5 and 1 m
            lift = min(lift, 0.5 * particle.depth)
            particle_ddepth -= lift/particle.fact        
        particle.last_update_res = time            

def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        print('lost Particle')
        particle.status = 10 #lost
        particle.delete()
        
def KeepInOcean(particle, fieldset, time):
    if particle.state == StatusCode.ErrorThroughSurface:
        particle.depth = 0.0
        particle.state = StatusCode.Success

def export(particle,fieldset,time):
    if particle.status==1:
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