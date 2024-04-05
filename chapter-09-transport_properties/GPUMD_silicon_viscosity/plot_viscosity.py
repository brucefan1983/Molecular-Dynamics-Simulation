import matplotlib.pyplot as plt
import numpy as np

viscosity=np.loadtxt('viscosity.out')

M=viscosity.shape[0]
t=viscosity[:,0] # ps

print('eta_L=',1000*np.mean(viscosity[125,10:13]),'+-',1000*np.std(viscosity[125,10:13])/np.sqrt(3))
print('eta_S=',1000*np.mean(viscosity[125,13:16]),'+-',1000*np.std(viscosity[125,13:16])/np.sqrt(3))

    
plt.figure(figsize=(8, 8))

plt.subplot(2, 1, 1)
plt.plot(t, 1000*viscosity[:,10], '-', linewidth=2,label='x')
plt.plot(t, 1000*viscosity[:,11], '-', linewidth=2,label='y')
plt.plot(t, 1000*viscosity[:,12], '-', linewidth=2,label='z')
plt.plot(t, 1000*np.mean(viscosity[:,10:13],axis=1), '--', linewidth=2,label='mean')
plt.xlabel('Correlation time (ps)', fontsize=15)
plt.ylabel('Longitudinal viscosity (mPa s)', fontsize=15)
plt.title('(a)', fontsize=15)
plt.legend(fontsize=15)


plt.subplot(2, 1, 2)
plt.plot(t, 1000*viscosity[:,13], '-', linewidth=2,label='xy')
plt.plot(t, 1000*viscosity[:,15], '-', linewidth=2,label='yz')
plt.plot(t, 1000*viscosity[:,14], '-', linewidth=2,label='zx')
plt.plot(t, 1000*np.mean(viscosity[:,13:16],axis=1), '--', linewidth=2,label='mean')
plt.xlabel('Correlation time (ps)', fontsize=15)
plt.ylabel('Shear viscosity (mPa s)', fontsize=15)
plt.title('(b)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-c9-viscosity.pdf')

