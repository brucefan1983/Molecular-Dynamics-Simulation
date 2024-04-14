import matplotlib.pyplot as plt
import numpy as np

dos=np.loadtxt('dos.out')
mvac=np.loadtxt('mvac.out')

print('normalization=',np.trapz(dos[:,1],dos[:,0]/2/np.pi)/13824)

plt.figure(figsize=(10, 5))
plt.subplot(1,2,1)
plt.plot(mvac[:,0],mvac[:,1], '-', linewidth=3,label='x')
plt.plot(mvac[:,0],mvac[:,2], '--', linewidth=2,label='y')
plt.plot(mvac[:,0],mvac[:,3], ':', linewidth=1,label='z')
plt.text(0.05,0.9,'(a)', fontsize=15)
plt.xlabel('Correlation time (ps)', fontsize=15)
plt.ylabel('Normalized MVAC', fontsize=15)
plt.legend()
plt.tight_layout()

plt.subplot(1,2,2)
plt.plot(dos[:,0]/2/np.pi,dos[:,1], '-', linewidth=3,label='x')
plt.plot(dos[:,0]/2/np.pi,dos[:,2], '--', linewidth=2,label='y')
plt.plot(dos[:,0]/2/np.pi,dos[:,3], ':', linewidth=1,label='z')
plt.text(1,2350,'(b)', fontsize=15)
plt.xlabel('$\omega/2\pi$ (THz)', fontsize=15)
plt.ylabel('VDOS (THz$^{-1}$)', fontsize=15)
plt.tight_layout()

plt.savefig('fig-c09-vdos.pdf')



