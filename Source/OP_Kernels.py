def Buoyancy(particle, fieldset, time):
    '''Stokes law settling velocity'''
    if particle.beached == 0: #Check particle is in the water column
        if particle.tau==0: #Check age particle is 0
            if ParcelsRandom.uniform(1e-5,1) < particle.fratio: 
                particle.ro = 300 #randomly assign a fraction of the particles a different density, in this case floating density (keep a fraction of MP afloat)          
            particle.diameter = ParcelsRandom.normalvariate(particle.diameter, particle.SDD) #Randomly assign a value of diameter inside the Bamfield mesocosm size dist
            particle.length = ParcelsRandom.normalvariate(particle.length, particle.SDL) #Same for length
            particle.tau = 4*fieldset.rorunoff[time, particle.depth, 49.57871, -123.020164] #Assign Fraser river outflow at deploting time to particle (Used to calculate MP/m3)
        d = particle.diameter # particle diameter
        l = particle.length # particle length
        #visc=1e-3 #average viscosity sea water 
        bath = fieldset.bathym[time, particle.depth, particle.lat, particle.lon]
        if  particle.depth > bath: #Check bathymetry to trap in the sediment if too deep.
            particle.beached = 3 #particle trapped in the sediment
        else:
            g = 9.8 #Gravity
            rhob=1080 #HBac density fixed
            Vcell=8.3e-19 #volume Hbac cell fixed
            t = fieldset.T[time, particle.depth, particle.lat, particle.lon] #Loading temperature from SSC
            ro = fieldset.R[time, particle.depth, particle.lat, particle.lon] #Loading density sw from SSC
            NN = particle.Nbac #Number of bacteria attached to MP
            th= (Vcell*NN)/(2.5*2*math.pi*(d/2)*l+(d/2)**2) #rough approximation of thickness biofilm
            rho=-1000-ro+(particle.ro*l*(d/2)**2+rhob*(2*(d/2)**2*th+l*th**2+2*th**3))/(l*(d/2)**2+2*(d/2)**2*th+l*th**2+2*th**3) #Total density considering biofilm
            #difference Density sea water and particle: LDPE (~920 kg/m3 ),PS (~150 kg/m3), PET (~1370 kg/m3). 
            d+=2*th #diameter considering biofilm
            l+=2*th #length considering biofilm
            visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296) #kinematic viscosity for Temp of SSC
            Ws= ((l/d)**-1.664)*0.079*((l**2)*g*(rho))/(visc) #sinking velocity considering density and dimensions change from biofouling
            particle.dz = Ws*particle.dt
        #if dz+particle.depth > 0:
        #    particle.depth += dz #Change particle depth according to WS
        #else:
        #    particle.depth = dz #Keep particle near surface. assuming reflexion in the boundary

def DeleteParticle(particle, fieldset, time):
    """Delete particle from OceanParcels simulation to avoid run failure
    """
    
    print(f'Particle {particle.id} lost !! [{particle.time}, {particle.depth}, {particle.lat}, {particle.lon}]')
    particle.delete()

def turb_mix(particle,fieldset,time):
    if particle.depth + 0.5 > bath: #Only calculate gradient of diffusion for particles deeper than 0.6 otherwise OP will check for particles outside the domain and remove it.
        Kzdz = 0
    else: 
        Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth+0.5, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]) #forward difference 
    dgrad = Kzdz*particle.dt
    if particle.depth+0.5*dgrad > 0.5:
        Kz = fieldset.vert_eddy_diff[time, particle.depth+0.5*dgrad, particle.lat, particle.lon] #Vertical diffusivity SSC  #
    else:
        Kz = fieldset.vert_eddy_diff[time, 0.5, particle.lat, particle.lon] #Vertical diffusivity SSC  #
    Rr = ParcelsRandom.uniform(-1, 1)
    Rrx = ParcelsRandom.uniform(-1, 1)
    Rry = ParcelsRandom.uniform(-1, 1)

    d_random = sqrt(3*2*Kz*particle.dt) * Rr
    dzs = dgrad + particle.dz
    particle.dz = 0
  
    Dlayer = 0.5*sqrt(Kz*particle.dt) #mixing layer dependant on Kz
    if dzs < 0:
        if d_random + dzs +particle.depth > bath: #randomly in the water column
            particle.depth = bath - Dlayer * ParcelsRandom.uniform(0, 1)
        elif d_random + particle.depth +dzs < 0:
                particle.depth = Dlayer * ParcelsRandom.uniform(0, 1) #Well mixed boundary layer
        else:
            particle.depth += d_random + dzs
    else:
        if d_random + dzs +particle.depth > bath: #randomly in the water column
            particle.beached = 3
        elif d_random +particle.depth +dzs < 0:
                particle.depth = Dlayer * ParcelsRandom.uniform(0, 1) #Well mixed boundary layer
        else:
            particle.depth += d_random + dzs
    
    
def Stokes_drift(particle, fieldset, time):
    '''Stokes drift'''  
    if particle.beached == 0:
        lat = particle.lat
        if lat > 48 and lat < 51: #Check that particle is inside WW3 data field
            deg2met = 111319.5
            latT = 0.6495 #cos(particle.lat*(math.pi/180))
            (us0, vs0, wl) = fieldset.stokes[time, particle.depth, particle.lat, particle.lon]
            k = (2*math.pi)/wl
            us = (us0*exp(-math.fabs(2*k*particle.depth)))/(deg2met*latT)
            vs = (vs0*exp(-math.fabs(2*k*particle.depth)))/deg2met
            particle.lon += us * particle.dt 
            particle.lat += vs * particle.dt
        


def AdvectionRK4_3D(particle, fieldset, time):
    if particle.beached == 0: #Check particle is in the water column
        if particle.dz > 0:
            if particle.dz+particle.depth > 0:
                particle.depth += particle.dz #Change particle depth according to WS
            else:
                particle.depth = math.fabs(particle.dz) #Keep particle near surface. assuming reflexion in the boundary

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

def turb_mix2(particle,fieldset,time):
    if particle.beached==0:
        if particle.depth + 0.5 > bath: #Only calculate gradient of diffusion for particles deeper than 0.6 otherwise OP will check for particles outside the domain and remove it.
            Kzdz = 0
        else: 
            Kzdz = 2*(fieldset.vert_eddy_diff[time, particle.depth+0.5, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]) #forward difference 
        dgrad = Kzdz*particle.dt
        if particle.depth+0.5*dgrad > 0.5:
            Kz = fieldset.vert_eddy_diff[time, particle.depth+0.5*dgrad, particle.lat, particle.lon] #Vertical diffusivity SSC  #
        else:
            Kz = fieldset.vert_eddy_diff[time, 0.5, particle.lat, particle.lon] #Vertical diffusivity SSC  #
        Rr = ParcelsRandom.uniform(-1, 1)
        Rrx = ParcelsRandom.uniform(-1, 1)
        Rry = ParcelsRandom.uniform(-1, 1)

        kh = particle.Kh
        d_random = sqrt(3*2*Kz*particle.dt) * Rr
        d_x = sqrt(3*2*kh*particle.dt) * Rrx
        d_y = sqrt(3*2*kh*particle.dt) * Rry 
        dzs = dgrad + particle.dz
        particle.dz = 0
        d_randomx = particle.lon + d_x/(deg2met*latT)
        d_randomy = particle.lat + d_y/deg2met
        (uch,vch) = fieldset.UV[time, 0.5, d_randomy, d_randomx]
        if particle.lat < 49.237 and particle.lon > -123.196 and particle.lat > 49.074:
            pass #Dont let particles beach inside the fraser river
        elif uch == 0 and vch == 0:
            particle.beached = 1
        else:
            particle.lat = d_randomy
            particle.lon = d_randomx

        Dlayer = 0.5*sqrt(Kz*particle.dt) #mixing layer dependant on Kz
        if dzs < 0:
            if d_random + dzs +particle.depth > bath: #randomly in the water colum01n
                particle.depth = bath - Dlayer * ParcelsRandom.uniform(0, 1)
            elif d_random + particle.depth +dzs < 0:
                particle.depth = Dlayer * ParcelsRandom.uniform(0, 1) #Well mixed boundary layer
            else:
                particle.depth += d_random + dzs
        else:
            if d_random + dzs +particle.depth > bath: #randomly in the water column
                particle.beached = 3
            elif d_random +particle.depth +dzs < 0:
                    particle.depth = Dlayer * ParcelsRandom.uniform(0, 1) #Well mixed boundary layer
            else:
                particle.depth += d_random + dzs


def Beaching(particle, fieldset, time):
    '''Beaching prob'''  
    if particle.beached == 0: #Check particle is in the water column       
        Tb = particle.Lb*86400 #timescale beaching in seconds
        x_offset = particle.Db/(deg2met*latT) #Checking distance x of possible beaching
        y_offset = particle.Db/deg2met #Checking distance y of possible beaching
        Pb = 1 - exp(-particle.dt/Tb)
        if particle.lat < 48.6 and particle.lon < -124.7 or particle.lat < 49.237 and particle.lon > -123.196 and particle.lat > 49.074:
            pass #Dont let particles beach inside the fraser river
        elif ParcelsRandom.uniform(0,1)<Pb:
            DWS1 = fieldset.U[time, 0.5, particle.lat+y_offset, particle.lon+x_offset] #particle.depth 0.5 check surface beach
            DWS2 = fieldset.U[time, 0.5, particle.lat-y_offset, particle.lon+x_offset]  
            DWS3 = fieldset.U[time, 0.5, particle.lat-y_offset, particle.lon-x_offset]
            DWS4 = fieldset.U[time, 0.5, particle.lat+y_offset, particle.lon-x_offset]
            if DWS1 == 0 or DWS2 == 0 or DWS3 == 0 or DWS4 == 0:
                particle.beached = 1

def Beaching2(particle, fieldset, time):
    '''Beaching prob'''  
    if particle.beached == 0: #Check particle is in the water column       
        Tb = particle.Lb*86400 #timescale beaching in seconds
        Pb = 1 - exp(-particle.dt/Tb)
        if particle.lat < 48.6 and particle.lon < -124.7 or particle.lat < 49.237 and particle.lon > -123.196 and particle.lat > 49.074:
            pass #Dont let particles beach inside the fraser river
        elif fieldset.coast_mask[time, 1, particle.lat, particle.lon] == 1:
            if ParcelsRandom.uniform(0,1)<Pb:
                print('beached')
                particle.beached = 1


def Unbeaching(particle, fieldset, time):
    '''Resuspension prob'''  
    if particle.beached == 1:        
        Ub = particle.Ub*86400  #timescale unbeaching in seconds
        Pr = 1 - exp(-particle.dt/Ub)
        if ParcelsRandom.uniform(0,1)<Pr:
            particle.beached = 0
        

def Biofilm(particle, fieldset, time):
    Nbac = particle.Nbac
    Nflag = particle.Nflag
    #Vcell=8.3e-13 #volume Hbac cell
    D = particle.diameter*1e2 
    L = particle.length*1e2
    th2= (Vcell*NN)/(2.5*2*math.pi*(D/2)*L+(D/2)**2) #rough approximation 
    D+= 2*th2
    L+= 2*th2
    ESRt = (((D**2)*3*L/2)**(1/3))/2
    Cb = 1.5e6 #aver Bacterial abundance in SoG /cm3 (S.W. Wilhelm et al., 2001)
    Cf = 1650 #Estimation Proportional to Hbacteria abundance (Gasol,1994)
    ###Cf = fieldset.microzooplankton[time, particle.depth, particle.lat, particle.lon]*4733.5 #conversion from mmolNm3 to cell/cm3
    Db = 1.83e-5 #Diffusion Bacteria (Kiorbe et al, 2003) cm2/s
    Df = 5.83-5 #Diffusion Het.Nanoflag (Kiorbe et al, 2003)
    detb = 2.83e-4 #detaching rate bacteria (Kiorbe et al, 2003)
    detf = 6.667e-5 #detaching rate Het.Nanoflag (Kiorbe et al, 2003)
    pp = fieldset.PPDIATNO3[time, particle.depth, particle.lat, particle.lon]+fieldset.PPPHYNO3[time, particle.depth, particle.lat, particle.lon]
    grb = pp*2.65 #conversion from PP to bacterial growth rate considering 20% of PP ends up as BP
    fcl = 8.33e-9 #clearence rate nanoflagelates (Kiorbe et al, 2003)
    Pf = (fcl/(1+fcl*3.22e-2*(Nbac))) #flagellate grazing coefficient
    af = Pf*1e-2
    Betab = Db/ESRt
    Betaf = Df/ESRt
    Ap = 2*math.pi*((D/2)*L+(D/2)**2) #Surface area of particle.
    particle.Nbac += (Betab*Cb*Ap + (grb - detb)*Nbac -Pf*Nbac*Nflag)*particle.dt
    particle.Nflag +=  (Betaf*Cf*Ap + af*Nbac*Nflag - detf*Nflag)*particle.dt
    if particle.Nbac < 0:
        particle.Nbac = 0
    if particle.Nflag < 0:
        particle.Nflag = 0
