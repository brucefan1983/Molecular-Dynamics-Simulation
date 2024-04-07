import matplotlib.pyplot as plt
import numpy as np

thermo=np.loadtxt('thermo.out')
E=np.sum(thermo[:,1:3],axis=1)
p=np.mean(thermo[:,3:6],axis=1)
V=np.prod(thermo[:,9:12],axis=1)
H=E+p*V/160.217663
M=thermo.shape[0]
num_T=7
K=M//7
H=E.reshape((num_T,K))

Cp=np.zeros((5,K//2))
for n in np.arange(1,6):
    a1=H[n-1,K//2:]
    a2=H[n+1,K//2:]
    Cp[n-1]=(a2-a1)/100/8.617343e-5/8000
Cp_ave=np.mean(Cp,axis=1)
Cp_err=np.std(Cp,axis=1)/np.sqrt(K)
    
plt.figure(figsize=(6, 6))

plt.errorbar(np.arange(2,num_T)*50,Cp_ave,yerr=Cp_err,linewidth=3)
plt.xlabel('Temperature (K)', fontsize=15)
plt.ylabel('Cp ($k_B$/atom)', fontsize=15)
#plt.ylim((0,3.2))
plt.tight_layout()
plt.savefig('fig-Cp.pdf')