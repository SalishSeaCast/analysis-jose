def Buoyancy(particle, fieldset, time):
    '''Stokes law settling velocity'''
    rhob=1080
    rhop=1350
    if particle.beached == 0:
        if particle.tau==0:
            if ParcelsRandom.uniform(0,1) < particle.fratio:
                particle.ro = 1025
            particle.diameter = ParcelsRandom.normalvariate(particle.diameter, particle.SDD)
            particle.length = ParcelsRandom.normalvariate(particle.length, particle.SDL)
            particle.tau = 4*fieldset.rorunoff[time, particle.depth, 49.57871, -123.020164] #fraser river outflow released every second
        d = particle.diameter # particle diameter
        l = particle.length # particle length
        visc=1e-3 #average viscosity sea water 
        z = particle.depth
        bath = 10*fieldset.mbathy[time, particle.depth, particle.lat, particle.lon]
        if  z > bath:
            particle.beached = 3 #particle trapped in the sediment
        else:
            g = 9.8
            #t = fieldset.T[time, particle.depth, particle.lat, particle.lon]
            ro = fieldset.R[time, particle.depth, particle.lat, particle.lon]
            NN = particle.Nbac
            Vcell=8.3e-19 #volume Hbac cell
            #ath = 2*math.pi
            #bth = math.pi*l + 4*math.pi*(d/2)
            #cth = 2*math.pi*(d/2)**2 + 2*math.pi*(d/2)*l
            #dth = -Vcell*NN
            #pth = -bth/(3*ath)
            #qth = pth**3 + (bth*cth-3*ath*dth)/(6*ath**2)   
            #rth = cth/(3*ath)
            #th=((qth + (qth**2 + (rth-pth**2)**3)**(1/2))**(1/3) + (qth - (qth**2 + (rth-pth**2)**3)**(1/2))**(1/3) + pth).real - 2.033e-20
            th= (Vcell*NN)/(2.5*2*math.pi*(d/2)*l+(d/2)**2) #rough approximation 
            particle.ro=(rhop*l*(d/2)**2+rhob*(2*(d/2)**2*th+l*th**2+2*th**3))/(l*(d/2)**2+2*(d/2)**2*th+l*th**2+2*th**3)
            dro = particle.ro-1000-ro  #difference Density sea water and particle: LDPE (~920 kg/m3 ),PS (~150 kg/m3), PET (~1370 kg/m3). 
            particle.diameter += 2*th
            particle.length += 2*th
            #visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296)
            Ws= ((l/d)**-1.664)*0.079*((l**2)*g*(dro))/(visc)
            dz = Ws*particle.dt
            particle.tau += 1
        if dz+z > 0:
            particle.depth += dz 
        else:
            particle.depth = 0.1
        
def turb_mix(particle,fieldset,time):
    wprime = 0
    Kz = fieldset.vert_eddy_diff[time, particle.depth, particle.lat, particle.lon]
    if particle.depth > 1.1:
        Kzdz = (fieldset.vert_eddy_diff[time, particle.depth+1, particle.lat, particle.lon]-fieldset.vert_eddy_diff[time, particle.depth-1, particle.lat, particle.lon])/2
    else:
        Kzdz = 0
    dW = ParcelsRandom.normalvariate(0, sqrt(particle.dt))
    wprime = Kzdz + (sqrt(2*Kz)*dW)/particle.dt
    dzp = wprime*particle.dt
    if dzp+particle.depth > 0:
            particle.depth += dzp 
    else:
        particle.depth = 0.1


def Stokes_drift(particle, fieldset, time):
    '''Stokes drift'''  
    if particle.beached == 0:
        lat = particle.lat
        if lat > 48 and lat < 51: #Check that particle is inside WW3 data field
            deg2met = 111319.5
            latT = 0.682
            z0 = particle.depth
            (us0, vs0, wl) = fieldset.stokes[time, particle.depth, particle.lat, particle.lon]
            k = (2*math.pi)/wl
            us = (us0*exp(-abs(2*k*z0)))/(deg2met*latT)
            vs = (vs0*exp(-abs(2*k*z0)))/deg2met
            particle.lon += us * particle.dt
            particle.lat += vs * particle.dt
        
def DeleteParticle(particle, fieldset, time):
    """Delete particle from OceanParcels simulation to avoid run failure
    """
    
    print(f'Particle {particle.id} lost !! [{particle.time}, {particle.depth}, {particle.lat}, {particle.lon}]')
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


def Beaching(particle, fieldset, time):
    '''Beaching prob'''  
    if particle.beached == 0:        
        Tb = particle.Lb*86400 
        x_offset = particle.Db/111319.5
        y_offset = particle.Db/(111319.5*0.682)     
        Pb = 1 - exp(-particle.dt/Tb)
        if particle.lat < 48.6 and particle.lon < -124.7 or particle.lat < 49.237 and particle.lon > -123.196 and particle.lat > 49.074:
            pass
        elif ParcelsRandom.uniform(0,1)<Pb:
            DWS1 = fieldset.U[time, 0.5, particle.lat+y_offset, particle.lon+x_offset] #particle.depth 0.5 check surface beach
            DWS2 = fieldset.U[time, 0.5, particle.lat-y_offset, particle.lon+x_offset]
            DWS3 = fieldset.U[time, 0.5, particle.lat-y_offset, particle.lon-x_offset]
            DWS4 = fieldset.U[time, 0.5, particle.lat+y_offset, particle.lon-x_offset]
            if DWS1 == 0 or DWS2 == 0 or DWS3 == 0 or DWS4 == 0:
                particle.beached = 1


def Unbeaching(particle, fieldset, time):
    '''Resuspension prob'''  
    if particle.beached == 1:        
        Ub = particle.Ub*86400  
        Pr = 1 - exp(-particle.dt/Ub)
        if ParcelsRandom.uniform(0,1)<Pr:
            particle.beached = 0
     
        

def Biofilm(particle, fieldset, time):
    Nbac = particle.Nbac
    Nflag = particle.Nflag
    D = particle.diameter*1e2
    L = particle.length*1e2
    ESRt = (((D**2)*3*L/2)**(1/3))/2
    Cb = 1.5e6 #Bacterial abundance in the strait of georgia /cm3
    Cf = fieldset.microzooplankton[time, particle.depth, particle.lat, particle.lon]*4733.5 #conversion from mmolNm3 to cell/cm3
    Db = 1.83e-5 #Diffusion Bacteria (Kiorbe et al, 2003) cm2/s
    Df = 5.83-5 #Diffusion Het.Nanoflag (Kiorbe et al, 2003)
    detb = 2.83e-4 #detaching rate bacteria (Kiorbe et al, 2003)
    detf = 6.667e-5 #detaching rate Het.Nanoflag (Kiorbe et al, 2003)
    chl = fieldset.diatoms[time, particle.depth, particle.lat, particle.lon]+fieldset.flagellates[time, particle.depth, particle.lat, particle.lon]
    grb = chl*1.333/86400
    fcl = 8.33e-9 #clearence rate nanoflagelates (Kiorbe et al, 2003)
    Pf = (fcl/(1+fcl*3.22e-2*(Nbac))) #flagellate grazing coefficient
    af = Pf*1e-2
    Betab = Db/ESRt
    Betaf = Df/ESRt
    Ap = 2*math.pi*((D/2)*L+(D/2)**2) #Surface area of particle.
    particle.Nbac += (Betab*Cb*Ap + (grb - detb)*Nbac -Pf*Nbac*Nflag)*particle.dt
    particle.Nflag +=  (Betaf*Cf*Ap + af*Nbac*Nflag - detf*Nflag)*particle.dt
    