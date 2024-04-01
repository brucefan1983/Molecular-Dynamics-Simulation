import matplotlib.pyplot as plt
import numpy as np

sdc=np.loadtxt('sdc.out')
msd=np.loadtxt('msd.out')

M=sdc.shape[0]
t=sdc[:,0] # ps

print('SDC from VAC = ', np.mean(sdc[-1,4:7]), ' +- ', np.std(sdc[-1,4:7])/np.sqrt(3))
print('SDC from MSD = ', np.mean(msd[-1,4:7]), ' +- ', np.std(msd[-1,4:7])/np.sqrt(3))

    
plt.figure(figsize=(8, 9))

plt.subplot(3, 1, 1)
plt.plot(t, sdc[:,4], '-', linewidth=3,label='x')
plt.plot(t, sdc[:,5], '-', linewidth=3,label='y')
plt.plot(t, sdc[:,6], '-', linewidth=3,label='z')
plt.plot(t, np.mean(sdc[:,4:7],axis=1), '--', linewidth=2,label='mean')
plt.xlabel('Correlation time (ps)', fontsize=15)
plt.ylabel('SDC from VAC (Ang$^2$/ps)', fontsize=15)
plt.title('(a)', fontsize=15)
plt.legend(fontsize=15)


plt.subplot(3, 1, 2)
plt.plot(t, msd[:,4], '-', linewidth=3,label='x')
plt.plot(t, msd[:,5], '-', linewidth=3,label='y')
plt.plot(t, msd[:,6], '-', linewidth=3,label='z')
plt.plot(t, np.mean(sdc[:,4:7],axis=1), '--', linewidth=2,label='mean')
plt.xlim((0,4))
plt.xlabel('Correlation time (ps)', fontsize=15)
plt.ylabel('SDC fromo MSD (Ang$^2$/ps)', fontsize=15)
plt.title('(b)', fontsize=15)
plt.legend(fontsize=15)

plt.subplot(3, 1, 3)
plt.plot(t, np.mean(sdc[:,4:7],axis=1), '-', linewidth=3,label='From VAC')
plt.plot(t, np.mean(msd[:,4:7],axis=1), '--', linewidth=3,label='From MSD')
plt.xlim((0,4))
plt.xlabel('Correlation time (ps)', fontsize=15)
plt.ylabel('SDC (Ang$^2$/ps)', fontsize=15)
plt.title('(c)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-c9-sdc.pdf')

