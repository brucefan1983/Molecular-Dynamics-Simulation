import matplotlib.pyplot as plt
import numpy as np

thermo=np.loadtxt('thermo.out')

M=thermo.shape[0]
t=np.arange(1,M+1)*0.002 # ps

    
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(t, thermo[:,3], '-', linewidth=2,label='x')
plt.plot(t, thermo[:,4], '-', linewidth=2,label='y')
plt.plot(t, thermo[:,5], '-', linewidth=2,label='z')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Pressure (GPa)', fontsize=15)
plt.title('(a)', fontsize=15)
plt.legend(fontsize=15)


plt.subplot(1, 2, 2)
plt.plot(t, thermo[:,9], '-', linewidth=2,label='x')
plt.plot(t, thermo[:,10], '-', linewidth=2,label='y')
plt.plot(t, thermo[:,11], '-', linewidth=2,label='z')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Length ($\mathrm{\AA}$)', fontsize=15)
plt.title('(b)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-c07-npt-examples-1.pdf')


