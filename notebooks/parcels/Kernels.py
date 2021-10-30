def Buoyancy(particle, fieldset, time):
    '''Stokes law settling velocity and critical sinking velocity'''
    Rp = 1370 #Density particle: LDPE (~920 kg/m3 ),PS (~150 kg/m3), PET (~1370 kg/m3).
    ESD=5e-5 #Size particle (ESD) equivalent spherical diameter
    #visc=1e-3 #average viscosity sea water
     
    z = particle.depth
    u = fieldset.U[time, particle.depth, particle.lat, particle.lon]
    if  u == 0:
        dz = 0
    else:
        
        a = -0.15
        b = 0.78
        k=4.5e-3
        pa=10.13
        ro0 = 1027
        t0 = 10
        s0 = 35
        g = 9.8
        t = fieldset.T[time, particle.depth, particle.lat, particle.lon]
        s = fieldset.S[time, particle.depth, particle.lat, particle.lon]
        p = ((g*z*ro0)/1e4)+pa
        ro = ro0 + a*(t-t0) + b*(s-s0) + k*p
        visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296)
        Ws= ((ESD**2)*g*(Rp-ro))/(18*visc)
        WWS =(11.68 + 0.1991*ESD*1e6 + 0.0004*(ESD*1e6)**2- 0.0993*(Rp-ro) + 0.0002*(Rp-ro)**2)/86400
        
    if Ws>WWS:
        dz = WWS*particle.dt
    else:
        dz = Ws*particle.dt
    if dz+z > 0:
        particle.depth += dz
    else:
        particle.depth = 0.5
        
def DeleteParticle(particle, fieldset, time):
    """Delete particle from OceanParcels simulation to avoid run failure
    """
    
    print(f'Particle {particle.id} lost !! [{particle.lon}, {particle.lat}, {particle.depth}, {particle.time}]')
    particle.delete()
    
