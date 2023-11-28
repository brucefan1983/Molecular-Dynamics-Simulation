import matplotlib.pyplot as plt
import numpy as np

thermo=np.loadtxt('thermo.out')
temperature=thermo[:,0]
M=thermo.shape[0]
D=100
N=M//D
num_data=np.zeros(N)
sigma=np.zeros(N)
ave=np.zeros(N)
for n in range(N):
    num_data[n]=D*(n+1)
    ave[n]=np.mean(temperature[0:D*(n+1)])
    sigma[n]=np.std(temperature[0:D*(n+1)])/np.sqrt(D*(n+1))
    
plt.figure(figsize=(8, 8))

plt.subplot(2, 1, 1)
plt.plot(np.arange(1,M+1)*0.2, temperature, '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Temperature (K)', fontsize=15)
plt.title('(a)', fontsize=15)

plt.subplot(2, 1, 2)
plt.plot(num_data,sigma,'-', linewidth=2)
plt.xlabel('$M$', fontsize=15)
plt.ylabel('$\sigma_M$ (K)', fontsize=15)
plt.title('(b)', fontsize=15)


plt.tight_layout()
plt.savefig('fig-chapter-12-sigma_M.pdf')


