def Advection(particle, fieldset, time):
    if particle.beached == 0: #Check particle is in the water column
        if particle.tau==0: #Check age particle is 0    
            if ParcelsRandom.uniform(1e-5,1) < particle.fratio: 
            #    pass#
                particle.surf = 1  
            '''Assign length and diameter to particles following input distribution''' 
            particle.diameter = ParcelsRandom.normalvariate(particle.diameter, particle.SDD) #Randomly assign a value of diameter inside the Bamfield mesocosm size dist
            particle.length = ParcelsRandom.normalvariate(particle.length, particle.SDL) #Same for length    
        particle.tau += particle.dt
        if particle.tau > particle.dtmax:
            particle.delete()
        if particle.surf == 1: #If particle is floating safe w sampling
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
        else:
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
            wa = (w1 + 2*w2 + 2*w3 + w4) / 6.
            particle.wa = wa
            particle.lon += (u1 + 2*u2 + 2*u3 + u4) / 6. * particle.dt
            particle.lat += (v1 + 2*v2 + 2*v3 + v4) / 6. * particle.dt
            particle.depth += wa * particle.dt

def Stokes_drift(particle, fieldset, time):
    """Stokes drift"""
    if particle.beached == 0:
        lat = particle.lat
        if lat > 48 and lat < 51: #Check that particle is inside WW3 data field
            deg2met = 111319.5
            latT = 0.6495 #cos(particle.lat*(math.pi/180))
            (us0, vs0, wl) = fieldset.stokes[time, 0, particle.lat, particle.lon]
            k = (2*math.pi)/wl
            us = (us0*exp(-math.fabs(2*k*particle.depth)))/(deg2met*latT)
            vs = (vs0*exp(-math.fabs(2*k*particle.depth)))/deg2met
            particle.lon += us * particle.dt 
            particle.lat += vs * particle.dt

def Buoyancy(particle, fieldset, time):
    """Stokes law calculating settling velocity"""
    if particle.beached == 0: #Check particle is in the water column   
        deps = max(particle.depth,0.51)
        bath = fieldset.Bathymetry[time, 0, particle.lat, particle.lon]  
        d = particle.diameter # particle diameter
        l = particle.length # particle length
        g = 9.8 #Gravity
        #? rhob=1080 #HBac density fixed
        #? Vcell=8.3e-19 #volume Hbac cell fixed bacilus avg
        t = fieldset.votemper[time, deps, particle.lat, particle.lon] #Loading temperature from SSC
        ro = fieldset.sigma_theta[time, deps, particle.lat, particle.lon] #Loading density sw from SSC
        NN = particle.Nbac #Number of bacteria attached to MP
        #? th= (Vcell*NN)/(5*math.pi*(d/2)*l+(d/2)**2) #rough approximation of thickness biofilm
        #? rho=-1000-ro+particle.ro*l*(d/2)**2 + rhob*(1-((l*(d/2)**2)/(l*(d/2)**2 + 2*th*(d/2)**2 + l*th**2 + 2*th**3))) #Total density considering biofilm (- sw density) )
        #? d+=2*th #diameter considering biofilm
        #? l+=2*th #length considering biofilm
        visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296) #kinematic viscosity for Temp of SSC
        Ws= ((l/d)**-1.664)*0.079*((l**2)*g*(particle.ro-1000-ro))/(visc)
        #? Ws= ((l/d)**-1.664)*0.079*((l**2)*g*(rho))/(visc) #sinking velocity considering density and dimensions change from biofouling
        particle.ws = Ws

# def Biofilm(particle, fieldset, time):
#     if particle.beached == 0:
    #     Nflag = particle.Nflag
    #     Vcell=8.3e-13 #volume Hbac cell
    #     D = particle.diameter*1e2 
    #     L = particle.length*1e2
    #     th2= (Vcell*NN)/(2.5*2*math.pi*(D/2)*L+(D/2)**2) #rough approximation 
    #     D+= 2*th2
    #     L+= 2*th2
    #     ESRt = (((D**2)*3*L/2)**(1/3))/2
    #     Cb = 1.5e6 #aver Bacterial abundance in SoG /cm3 (S.W. Wilhelm et al., 2001)
    #     Cf = 1650 #Estimation Proportional to Hbacteria abundance (Gasol,1994) close to 1640 coastal surface normal average abundance. Fukami 1996 
    #     ###Cf = fieldset.microzooplankton[time, deps, particle.lat, particle.lon]*4733.5 #conversion from mmolNm3 to cell/cm3
    #     Db = 2.33e-5 #Diffusion Bacteria (Kiorbe et al, 2003) cm2/s
    #     Df = 9.8e-5 #Diffusion Het.Nanoflag (Kiorbe et al, 2003)
    #     detb = 2.83e-4 #detaching rate bacteria (Kiorbe et al, 2003)
    #     detf = 6.667e-5 #detaching rate Het.Nanoflag (Kiorbe et al, 2003)
    #     pp = fieldset.PPDIATNO3[time, deps, particle.lat, particle.lon]+fieldset.PPPHYNO3[time, deps, particle.lat, particle.lon]
    #     grb = pp*2.65 #conversion from PP to bacterial growth rate considering 20% of PP ends up as BP
    #     fcl = 8.33e-9 #clearence rate nanoflagelates (Kiorbe et al, 2003)
    #     Pf = (fcl/(1+fcl*3.22e-2*(NN))) #flagellate grazing coefficient
    #     af = Pf*1e-2
    #     Betab = Db/(ESRt*100)
    #     Betaf = Df/(ESRt*100)
    #     Ap = 2*math.pi*((D/2)*L+(D/2)**2) #Surface area of particle.
    #     particle.Nbac += (Betab*Cb*Ap + (grb - detb)*NN -Pf*NN*Nflag)*particle.dt
    #     particle.Nflag +=  (Betaf*Cf*Ap + af*NN*Nflag - detf*Nflag)*particle.dt
    #     if particle.Nbac < 0:
    #         particle.Nbac = 0
    #     if particle.Nflag < 0:
    #         particle.Nflag = 0            

def turb_mix(particle,fieldset,time):
    """Vertical mixing"""
    if particle.beached==0:
        #Vertical mixing
        if particle.depth + 0.5 > bath: #Only calculate gradient of diffusion for particles deeper than 0.6 otherwise OP will check for particles outside the domain and remove it.
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth-0.5, particle.lat, particle.lon]) #backwards difference 
        else: 
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth+0.5, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]) #forward difference 
        dgrad = Kzdz*particle.dt
        if particle.depth+0.5*dgrad > 0:
            Kz = fieldset.vert_eddy_diff[time, particle.depth+0.5*dgrad, particle.lat, particle.lon] #Vertical diffusivity SSC  
        else:
            Kz = fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon] 
        if particle.depth+0.5*dgrad > bath:  
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
        Sbh = fieldset.vosaline[time, deps, d_randomy, d_randomx] #Check if particles reach coast at surface (Salinity = 0)
    
def Displacement(particle,fieldset,time):
    ''''Apply movement calculated by other kernels'''
    if particle.beached==0 and particle.surf == 0:
        #Apply turbulent mixing.
        if dzs + particle.depth > bath: #randomly in the water column
            particle.depth = bath - Dlayer * ParcelsRandom.uniform(0, 1) #Well mixed boundary layer
        elif particle.depth + dzs < 0:
            #if ParcelsRandom.uniform(1e-5,1) < particle.fratio:
            #    particle.surf = 1
            #    print('surfaced')
            #else:
            particle.depth = Dlayer * ParcelsRandom.uniform(0, 1) #Well mixed boundary layer
        else:
            particle.depth += dzs #apply mixing
        #Apply horizontal mixing (beaching for particles pushed through coast) 
        if particle.lat < 49.237 and particle.lon > -123.196 and particle.lat > 49.074:
            pass #Dont let particles beach inside the Fraser river
        elif Sbh == 0:
            particle.beached = 1
        else:
            particle.lat = d_randomy
            particle.lon = d_randomx
        #Apply buoyancy
        if particle.depth + Ws*particle.dt > bath:
            particle.beached = 3 #Trap particle in sediment (sticky bottom)
        elif particle.depth + Ws*particle.dt < 0: 
            particle.depth = 0
        else:
            particle.depth += Ws*particle.dt


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

def DeleteParticle(particle, fieldset, time):
    """Delete particle from OceanParcels simulation to avoid run failure"""
    #particle.beached = 4
    print(f'Particle {particle.id} lost !! [{particle.time}, {particle.depth}, {particle.lat}, {particle.lon}]')
    particle.delete()


# def Stokes_driftRK3(particle, fieldset, time):
#     """Stokes drift solved with 3rd order RungeKutta"""
#     if particle.beached == 0:
#         lat = particle.lat
#         if lat > 48 and lat < 51: #Check that particle is inside WW3 data field
#             deg2met = 111319.5
#             latT = 0.6495 #cos(particle.lat*(math.pi/180))
#             (us0, vs0, wl) = fieldset.stokes[time, particle.depth, particle.lat, particle.lon]
#             k = (2*math.pi)/wl
#             us1 = (us0*exp(-math.fabs(2*k*particle.depth)))/(deg2met*latT)
#             vs1 = (vs0*exp(-math.fabs(2*k*particle.depth)))/deg2met
#             lon1s = particle.lon + us1 * (1/3) *particle.dt 
#             lat1s = particle.lat + vs1 * (1/3) * particle.dt
#             (us0, vs0, wl) = fieldset.stokes[time + (1/3) * particle.dt, particle.depth, lat1s, lon1s]
#             k = (2*math.pi)/wl
#             us2 = (us0*exp(-math.fabs(2*k*particle.depth)))/(deg2met*latT)
#             vs2 = (vs0*exp(-math.fabs(2*k*particle.depth)))/deg2met
#             lon2s = particle.lon + us2 * .5 *particle.dt 
#             lat2s = particle.lat + vs2 * .5 * particle.dt
#             (us0, vs0, wl) = fieldset.stokes[time + .5 * particle.dt, particle.depth, lat2s, lon2s]
#             k = (2*math.pi)/wl
#             us3 = (us0*exp(-math.fabs(2*k*particle.depth)))/(deg2met*latT)
#             vs3 = (vs0*exp(-math.fabs(2*k*particle.depth)))/deg2met
#             particle.lon +=  us3 * .5 *particle.dt 
#             particle.lat +=  vs3 * .5 * particle.dt