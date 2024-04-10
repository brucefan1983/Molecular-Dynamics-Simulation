import matplotlib.pyplot as plt
import numpy as np

M=1000
data=np.loadtxt('compute.out')
dEdt1=data[-1,-1]/M; #eV/ps
dEdt2=-data[-1,-2]/M; #eV/ps
A=6.57204**2 # nm^2 
dT=np.mean(data[500:,1])-np.mean(data[500:,5]); # K
G=np.array((1.602177e+2*dEdt1/A/dT,1.602177e+2*dEdt2/A/dT))

print('G=',np.mean(G),'+-',np.std(G)/np.sqrt(2))

plt.figure(figsize=(8, 8))
plt.subplot(2,1,1)
plt.plot(np.arange(1,6),np.mean(data[500:,1:6],axis=0), '-o', linewidth=2)
plt.title('(a)')
plt.xlabel('Group index', fontsize=15)
plt.ylabel('Temperature (K)', fontsize=15)
plt.tight_layout()

plt.subplot(2,1,2)
plt.plot(np.arange(M),data[:,-2], '-', linewidth=2,label='source')
plt.plot(np.arange(M),data[:,-1], '--', linewidth=2,label='sink')
plt.title('(b)')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Energy (eV)', fontsize=15)
plt.legend()
plt.tight_layout()

plt.savefig('fig-c9-kappa-nemd.pdf')