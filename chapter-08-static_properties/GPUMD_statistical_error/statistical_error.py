import matplotlib.pyplot as plt
import numpy as np

thermo=np.loadtxt('thermo.out')
temperature=thermo[:,0]
pressure=np.mean(thermo[:,3:6],axis=1)
M=thermo.shape[0]
D=100
N=M//D
num_data=np.zeros(N)
T_sigma=np.zeros(N)
T_ave=np.zeros(N)
P_sigma=np.zeros(N)
P_ave=np.zeros(N)
for n in range(N):
    num_data[n]=D*(n+1)
    T_ave[n]=np.mean(temperature[0:D*(n+1)])
    T_sigma[n]=np.std(temperature[0:D*(n+1)])/np.sqrt(D*(n+1))
    P_ave[n]=np.mean(pressure[0:D*(n+1)])
    P_sigma[n]=np.std(pressure[0:D*(n+1)])/np.sqrt(D*(n+1))
    
plt.figure(figsize=(8, 8))

plt.subplot(2, 1, 1)
plt.plot(np.arange(1,M+1)*0.2, temperature, '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Temperature (K)', fontsize=15)
plt.title('(a)', fontsize=15)

plt.subplot(2, 1, 2)
plt.plot(num_data,T_sigma,'-', linewidth=2)
plt.xlabel('$M$', fontsize=15)
plt.ylabel('$\sigma_M$ (K)', fontsize=15)
plt.title('(b)', fontsize=15)

plt.tight_layout()
plt.savefig('fig-chapter-12-sigma_T.pdf')

plt.figure(figsize=(8, 8))

plt.subplot(2, 1, 1)
plt.plot(np.arange(1,M+1)*0.2, pressure, '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Pressure (GPa)', fontsize=15)
plt.title('(a)', fontsize=15)

plt.subplot(2, 1, 2)
plt.plot(num_data,P_sigma,'-', linewidth=2)
plt.xlabel('$M$', fontsize=15)
plt.ylabel('$\sigma_M$ (GPa)', fontsize=15)
plt.title('(b)', fontsize=15)

plt.tight_layout()
plt.savefig('fig-chapter-12-sigma_P.pdf')


