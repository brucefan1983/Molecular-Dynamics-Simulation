import matplotlib.pyplot as plt
import numpy as np

thermo_ber=np.loadtxt('ber/thermo.out')
thermo_scr=np.loadtxt('scr/thermo.out')
thermo_mttk=np.loadtxt('mttk/thermo.out')

M=thermo_ber.shape[0]
t=np.arange(1,M+1)*0.001 # ps

    
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(t, np.mean(thermo_mttk[:,3:5],axis=1), '-', linewidth=2,label='MTTK')
plt.plot(t, np.mean(thermo_scr[:,3:5],axis=1), '-', linewidth=2,label='SCR')
plt.plot(t, np.mean(thermo_ber[:,3:5],axis=1), '-', linewidth=2,label='Berendsen')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Pressure (GPa)', fontsize=15)
plt.title('(a)', fontsize=15)
plt.legend(fontsize=15)


plt.subplot(1, 2, 2)
plt.plot(t, np.mean(thermo_mttk[:,9:12],axis=1), '-', linewidth=2,label='MTTK')
plt.plot(t, np.mean(thermo_scr[:,9:12],axis=1), '-', linewidth=2,label='SCR')
plt.plot(t, np.mean(thermo_ber[:,9:12],axis=1), '-', linewidth=2,label='Berendsen')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Length ($\mathrm{\AA}$)', fontsize=15)
plt.title('(b)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-c07-barostat-compare.pdf')


