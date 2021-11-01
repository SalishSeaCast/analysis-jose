import numpy as np
import matplotlib.pyplot as plt

g = 9.8


Rw = 1027
vk = 1e-3/Rw

Dp = np.linspace(5e-5,5e-3,20)
Rp = np.linspace(800,1500,20)
Dc = []
Wb1 = np.zeros([20,20])
Wb2 = np.zeros([20,20])

for i,Ro in enumerate(Rp):
    Dc.append((9.52*(vk**0.666))/((g**0.333)*(1-((Ro/Rw)**0.333))))
    for j,D in enumerate(Dp):
        Wb1[i,j] = ((g*(D**2)*(1-(Ro/Rw)))/(18*vk))
        Wb2[i,j] = (8/3)*g*D*np.sign(1-(Ro/Rw))*np.sqrt(np.abs(1-(Ro/Rw)))

#plt.plot(Rp,Dc)
plt.plot(Dp,Wb1[10,:])
plt.plot(Dp,Wb2[10,:])
plt.show()