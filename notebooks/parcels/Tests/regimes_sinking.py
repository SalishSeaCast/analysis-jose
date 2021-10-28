import numpy as np
from matplotlib import pyplot as plt

n=100
Rp = np.linspace(1000,1500,n) #Density particle: LDPE (~920 kg/m3 ),PS (~150 kg/m3), PET (~1370 kg/m3).
ESD=5e-5 #Size particle (ESD) equivalent spherical diameter
#visc=1e-3 #average viscosity sea water



ro = 1027
g = 9.8
t=18

visc = 4.2844e-5 + 1/(0.157*((t + 64.993)**2)-91.296)

Ws= ((ESD**2)*g*(Rp/ro-1))/(18*visc)*86400
WWS =(11.68 + 0.1991*ESD*1e6 + 0.0004*(ESD*1e6)**2- 0.0993*(Rp-ro) + 0.0002*(Rp-ro)**2)

plt.plot(Rp,Ws,'b')
plt.plot(Rp,WWS,'r')
plt.plot(Rp,np.zeros(n),'k--')
plt.axvline(x=ro,linewidth=1, color='b')
plt.axvline(x=920,linewidth=1, color='g')
plt.axvline(x=150,linewidth=1, color='g')
plt.axvline(x=1370,linewidth=1, color='g')
plt.show()
