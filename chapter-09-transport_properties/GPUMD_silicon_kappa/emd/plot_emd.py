import matplotlib.pyplot as plt
import numpy as np

hac=np.loadtxt('hac.out')

M=5000
Ns=hac.shape[0]//M
t=hac[0:M,0] # ps
kappa=np.sum(hac[:,6:],axis=1)/3
kappa=kappa.reshape((Ns,M))

print('Ns=',Ns)
print('kappa=',np.mean(kappa[:,-1]),'+-',np.std(kappa[:,-1])/np.sqrt(Ns))

    
plt.figure(figsize=(8, 5))
plt.plot(t, kappa.transpose(), '--',color=(0.7,0.7,0.7), linewidth=1)
plt.plot(t, np.mean(kappa,axis=0), '-', linewidth=3)
plt.xlabel('Correlation time (ps)', fontsize=20)
plt.ylabel('$\kappa$ (W m$^{-1}$ K$^{-1}$)', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
plt.savefig('fig-c09-kappa-emd.pdf')


