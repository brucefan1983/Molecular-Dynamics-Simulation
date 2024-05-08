import matplotlib.pyplot as plt
import numpy as np

kappa=np.loadtxt('kappa.out')

M=10000
Ns=kappa.shape[0]//M
t=np.arange(M) # ps
kappa=np.sum(kappa[:,0:2],axis=1)
kappa=kappa.reshape((Ns,M))
kappa_cum=np.zeros((Ns,M))
for ns in range(Ns):
    kappa_cum[ns,:] = np.cumsum(kappa[ns,:])/np.arange(1,M+1)

print('Ns=',Ns)
print('kappa=',np.mean(kappa_cum[:,-1]),'+-',np.std(kappa_cum[:,-1])/np.sqrt(Ns))

    
plt.figure(figsize=(12, 6))
plt.subplot(1,2,1)
plt.plot(t, kappa[0,:], '-',color=(0.7,0.7,0.7), linewidth=1)
plt.plot(t, kappa_cum[0,:], '--', linewidth=1)
plt.title('(a)', fontsize=25)
plt.xlabel('Time (ps)', fontsize=25)
plt.ylabel('$\kappa$ (W m$^{-1}$ K$^{-1}$)', fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.subplot(1,2,2)

plt.plot(t, kappa_cum.transpose(), '--', linewidth=1)
plt.plot(t, np.mean(kappa_cum,axis=0), '-', linewidth=3)
plt.title('(b)', fontsize=25)
plt.xlabel('Time (ps)', fontsize=25)
plt.ylabel('$\kappa$ (W m$^{-1}$ K$^{-1}$)', fontsize=25)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()

plt.savefig('fig-c09-kappa-hnemd.pdf')