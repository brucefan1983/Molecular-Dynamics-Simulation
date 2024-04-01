import matplotlib.pyplot as plt
import numpy as np

msd=np.loadtxt('msd.out')

M=msd.shape[0]
t=msd[:,0] # ps

    
plt.figure(figsize=(8, 6))

plt.plot(t, msd[:,1], '-', linewidth=2,label='x')
plt.plot(t, msd[:,2], '-', linewidth=2,label='y')
plt.plot(t, msd[:,3], '-', linewidth=2,label='z')
plt.plot(t, np.mean(msd[:,1:4],axis=1), '--', linewidth=2,label='mean')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('MSD (Ang$^2$)', fontsize=15)
plt.legend(fontsize=15)


plt.tight_layout()
plt.savefig('fig-c9-msd.pdf')

