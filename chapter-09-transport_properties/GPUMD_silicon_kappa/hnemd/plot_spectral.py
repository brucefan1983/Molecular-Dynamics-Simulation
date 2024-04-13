import matplotlib.pyplot as plt
import numpy as np

shc=np.loadtxt('shc.out')
Ns=4
M=shc.shape[0]//Ns
shc_ave=np.zeros((M,3))
for ns in range(Ns):
    shc_ave=shc_ave+shc[ns*M:(ns+1)*M,:]
shc_ave=shc_ave/4

Nc=250
Nt=2*Nc-1
a=5.4767 # length of the group in units of Angstrom
V=a**3*12*12*12
t=shc_ave[0:Nt,0]
K=np.sum(shc_ave[0:Nt,1:3],axis=1) # eV A/ps
omega=shc_ave[Nt:,0] 
nu=omega/2/np.pi
T=300
Fe=1.5e-5
k=np.sum(shc_ave[Nt:,1:3],axis=1)*1.602177e+3/(V*T*Fe)

print('k=',np.sum(k)*(omega[1]-omega[0])/2/np.pi)

plt.figure(figsize=(8, 8))
plt.subplot(2,1,1)
plt.plot(t,K, '-', linewidth=2)
plt.title('(a)')
plt.xlabel('Correlation time (ps)', fontsize=15)
plt.ylabel('K (eV Ang ps$^{-1}$)', fontsize=15)
plt.tight_layout()

plt.subplot(2,1,2)
plt.plot(omega/2/np.pi,k, '-', linewidth=2)
plt.title('(b)')
plt.xlabel('$\omega/2\pi$ (THz)', fontsize=15)
plt.ylabel('G (GW m$^{-2}$ K$^{-1}$ THz$^{-1}$)', fontsize=15)
plt.tight_layout()

plt.savefig('fig-c9-kappa-k-omega.pdf')


np.savetxt('nu.txt',nu)
np.savetxt('k.txt',k)


