import matplotlib.pyplot as plt
import numpy as np

thermo_ber=np.loadtxt('ber/thermo.out')
thermo_scr=np.loadtxt('scr/thermo.out')
thermo_mttk=np.loadtxt('mttk/thermo.out')

M=thermo_ber.shape[0]
t=np.arange(1,M+1)*0.001 # ps

plt.figure()
plt.plot(t, np.mean(thermo_mttk[:,3:5],axis=1), '-', linewidth=2,label='MTTK')
plt.plot(t, np.mean(thermo_scr[:,3:5],axis=1), '-', linewidth=2,label='SCR')
plt.plot(t, np.mean(thermo_ber[:,3:5],axis=1), '-', linewidth=2,label='Berendsen')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Pressure (GPa)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-barostat-compare.pdf')

plt.figure()
plt.plot(t, np.mean(thermo_mttk[:,9:12],axis=1), '-', linewidth=2,label='MTTK')
plt.plot(t, np.mean(thermo_scr[:,9:12],axis=1), '-', linewidth=2,label='SCR')
plt.plot(t, np.mean(thermo_ber[:,9:12],axis=1), '-', linewidth=2,label='Berendsen')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Length (A)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-barostat-compare-length.pdf')



plt.figure()
plt.hist(np.mean(thermo_scr[500:,3:5],axis=1),label='SCR',bins=51)
plt.hist(np.mean(thermo_mttk[500:,3:5],axis=1),label='MTTK',bins=51)

#plt.hist(np.mean(thermo_ber[:,3:5],axis=1),label='Berendsen',bins=51)
plt.ylabel('Count', fontsize=15)
plt.xlabel('Pressure (GPa)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-barostat-compare2.pdf')


