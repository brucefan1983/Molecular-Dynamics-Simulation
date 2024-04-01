import matplotlib.pyplot as plt
import numpy as np

sdc_nve=np.loadtxt('nve/sdc.out')
sdc_bdp=np.loadtxt('bdp/sdc.out')
sdc_lan=np.loadtxt('lan/sdc.out')

M=sdc_nve.shape[0]
t=sdc_nve[:,0] # ps

print('SDC from NVE = ', np.mean(sdc_nve[-1,4:7]), ' +- ', np.std(sdc_nve[-1,4:7])/np.sqrt(3))
print('SDC from BDP = ', np.mean(sdc_bdp[-1,4:7]), ' +- ', np.std(sdc_bdp[-1,4:7])/np.sqrt(3))
print('SDC from LAN = ', np.mean(sdc_lan[-1,4:7]), ' +- ', np.std(sdc_lan[-1,4:7])/np.sqrt(3))
    
plt.figure(figsize=(8, 6))

plt.plot(t, np.mean(sdc_nve[:,4:7],axis=1), '-', linewidth=3,label='NVE')
plt.plot(t, np.mean(sdc_bdp[:,4:7],axis=1), '-', linewidth=3,label='NVT-BDP')
plt.plot(t, np.mean(sdc_lan[:,4:7],axis=1), '-', linewidth=3,label='NVT-Langevin')
plt.xlim((0,4))
plt.xlabel('Correlation time (ps)', fontsize=15)
plt.ylabel('SDC (Ang$^2$/ps)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-c9-diffusion-thermostat.pdf')

