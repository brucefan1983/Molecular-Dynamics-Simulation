import matplotlib.pyplot as plt
import numpy as np

thermo_ber=np.loadtxt('ber/thermo.out')
thermo_bdp=np.loadtxt('bdp/thermo.out')
thermo_nhc=np.loadtxt('nhc/thermo.out')
thermo_lan=np.loadtxt('lan/thermo.out')

M=thermo_ber.shape[0]
t=np.arange(1,M+1)*0.001 # ps

plt.figure()
plt.plot(t, thermo_ber[:,0], '-', linewidth=2,label='Berendsen')
plt.plot(t, thermo_bdp[:,0], '-', linewidth=2,label='BDP')
plt.plot(t, thermo_nhc[:,0], '-', linewidth=2,label='NHC')
plt.plot(t, thermo_lan[:,0], '-', linewidth=2,label='Langevin')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Temperature (K)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-chapter-10-compare-thermostats.pdf')


