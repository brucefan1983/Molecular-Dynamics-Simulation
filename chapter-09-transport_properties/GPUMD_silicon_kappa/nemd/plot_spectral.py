import matplotlib.pyplot as plt
import numpy as np

data=np.loadtxt('compute.out')
dT=np.mean(data[500:,1])-np.mean(data[500:,5]); # K

Nc=250
Nt=2*Nc-1
a=5.4767 # length of the group in units of Angstrom
V=a**3*12*12
shc=np.loadtxt('shc.out')
t=shc[0:Nt,0]
K=np.sum(shc[0:Nt,1:3],axis=1) # eV A/ps
omega=shc[Nt:,0] 
nu=omega/2/np.pi
G=np.sum(shc[Nt:,1:3],axis=1)*1.602177e+4/V/dT

print('G=',np.sum(G)*(omega[1]-omega[0])/2/np.pi)

plt.figure(figsize=(8, 8))
plt.subplot(2,1,1)
plt.plot(t,K, '-', linewidth=2)
plt.title('(a)')
plt.xlabel('Correlation time (ps)', fontsize=15)
plt.ylabel('K (eV Ang ps$^{-1}$)', fontsize=15)
plt.tight_layout()

plt.subplot(2,1,2)
plt.plot(omega/2/np.pi,G, '-', linewidth=2)
plt.title('(b)')
plt.xlabel('$\omega/2\pi$ (THz)', fontsize=15)
plt.ylabel('G (GW m$^{-2}$ K$^{-1}$ THz$^{-1}$)', fontsize=15)
plt.tight_layout()

plt.savefig('fig-c9-kappa-G-omega.pdf')

np.savetxt('nu.txt',nu)
np.savetxt('G.txt',G)
