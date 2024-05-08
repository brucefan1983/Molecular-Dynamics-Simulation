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

plt.figure(figsize=(12, 6))
plt.subplot(1,2,1)
plt.plot(t,K, '-', linewidth=3)
plt.title('(a)', fontsize=25)
plt.xlabel('Correlation time (ps)', fontsize=25)
plt.ylabel('K (eV $\mathrm{\AA}$ ps$^{-1}$)', fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(omega/2/np.pi,G, '-', linewidth=3)
plt.title('(b)', fontsize=25)
plt.xlabel('$\omega/2\pi$ (THz)', fontsize=25)
plt.ylabel('G (GW m$^{-2}$ K$^{-1}$ THz$^{-1}$)', fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.savefig('fig-c09-kappa-G-omega.pdf')

np.savetxt('nu.txt',nu)
np.savetxt('G.txt',G)
