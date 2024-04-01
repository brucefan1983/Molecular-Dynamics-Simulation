import matplotlib.pyplot as plt
import numpy as np

sdc=np.loadtxt('nve/sdc.out')

M=sdc.shape[0]
t=sdc[:,0] # ps

    
plt.figure(figsize=(8, 8))

plt.subplot(2, 1, 1)
plt.plot(t, sdc[:,1], '-', linewidth=2,label='x')
plt.plot(t, sdc[:,2], '-', linewidth=2,label='y')
plt.plot(t, sdc[:,3], '-', linewidth=2,label='z')
plt.plot(t, np.mean(sdc[:,1:4],axis=1), '--', linewidth=2,label='mean')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('VAC (Ang$^2$/ps$^2$)', fontsize=15)
plt.legend(fontsize=15)


plt.subplot(2, 1, 2)
plt.plot(t, sdc[:,1], '-', linewidth=2,label='x')
plt.plot(t, sdc[:,2], '-', linewidth=2,label='y')
plt.plot(t, sdc[:,3], '-', linewidth=2,label='z')
plt.plot(t, np.mean(sdc[:,1:4],axis=1), '--', linewidth=2,label='mean')
plt.xlim((0,0.4))
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('VAC (Ang$^2$/ps$^2$)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-c9-vac.pdf')

