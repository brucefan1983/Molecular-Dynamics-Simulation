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
plt.plot(omega/2/np.pi,k, '-', linewidth=3)
plt.title('(b)', fontsize=25)
plt.xlabel('$\omega/2\pi$ (THz)', fontsize=25)
plt.ylabel('G (GW m$^{-2}$ K$^{-1}$ THz$^{-1}$)', fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.savefig('fig-c09-kappa-k-omega.pdf')


np.savetxt('nu.txt',nu)
np.savetxt('k.txt',k)


