def Buoyancy(particle, fieldset, time):
    '''Stokes law settling velocity'''
    Rp = particle.ro #Density particle: LDPE (~920 kg/m3 ),PS (~150 kg/m3), PET (~1370 kg/m3).
    ESD= particle.size #Size particle (ESD) equivalent spherical diameter
    visc=1e-3 #average viscosity sea water
     
    z = particle.depth
    u = fieldset.U[time, particle.depth, particle.lat, particle.lon]
    if  u == 0:
        dz = 0
    else:
        g = 9.8
        #t = fieldset.T[time, particle.depth, particle.lat, particle.lon]
        ro = fieldset.R[time, particle.depth, particle.lat, particle.lon]
        #visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296)
        Ws= ((ESD**2)*g*(Rp-ro))/(18*visc)
        WWS =(11.68 + 0.1991*ESD*1e6 + 0.0004*(ESD*1e6)**2- 0.0993*(Rp-ro) + 0.0002*(Rp-ro)**2)/86400
        if Ws > WWS:
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
    
