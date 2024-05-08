import matplotlib.pyplot as plt
import numpy as np

thermo=np.loadtxt('thermo.out')

M=thermo.shape[0]
t=np.arange(1,M+1)*0.002 # ps

    
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(t, thermo[:,3], '-', linewidth=2,label='xx')
plt.plot(t, thermo[:,4], '-', linewidth=2,label='yy')
plt.plot(t, thermo[:,5], '-', linewidth=2,label='zz')
plt.plot(t, thermo[:,6], '-', linewidth=2,label='yz')
plt.plot(t, thermo[:,7], '-', linewidth=2,label='zx')
plt.plot(t, thermo[:,8], '--', linewidth=2,label='xy')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Pressure (GPa)', fontsize=15)
plt.title('(a)', fontsize=15)
plt.legend(fontsize=15)


plt.subplot(1, 2, 2)
plt.plot(t, thermo[:,10], '-', linewidth=2,label='a_y')
plt.plot(t, thermo[:,12], '--', linewidth=2,label='b_x')
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Length ($\mathrm{\AA}$)', fontsize=15)
plt.title('(b)', fontsize=15)
plt.legend(fontsize=15)

plt.tight_layout()
plt.savefig('fig-c07-npt-examples-2.pdf')


