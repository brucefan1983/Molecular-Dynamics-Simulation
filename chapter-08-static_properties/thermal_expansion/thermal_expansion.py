import matplotlib.pyplot as plt
import numpy as np

thermo=np.loadtxt('thermo.out')
temperature=thermo[:,0]
pressure=np.mean(thermo[:,3:6],axis=1)
box_length=np.mean(thermo[:,9:12],axis=1)
M=thermo.shape[0]
num_T=7
K=M//7
t=np.arange(1,M+1)*0.2
a=box_length.reshape((num_T,K))/10
a_ave=np.mean(a[:,K//2:],axis=1)

plt.figure(figsize=(8, 8))

plt.subplot(2, 2, 1)
plt.plot(t,temperature, '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Temperature (K)', fontsize=15)
plt.title('(a)', fontsize=15)

plt.subplot(2, 2, 2)
plt.plot(t,pressure, '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Pressure (GPa)', fontsize=15)
plt.title('(b)', fontsize=15)

plt.subplot(2, 2, 3)
plt.plot(t,box_length, '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Box length (A)', fontsize=15)
plt.title('(c)', fontsize=15)

plt.subplot(2, 2, 4)
plt.plot(np.arange(1,num_T+1)*50,a_ave, 'o', linewidth=2)
plt.xlabel('Temperature (K)', fontsize=15)
plt.ylabel('a (A)', fontsize=15)
plt.title('(d)', fontsize=15)

plt.tight_layout()
plt.savefig('fig-chapter-12-expansion1.pdf')

expansion=np.zeros((5,K//2))
for n in np.arange(1,6):
    a1=a[n-1,K//2:]
    a2=a[n+1,K//2:]
    expansion[n-1]=(a2-a1)*2/(a2+a1)/100
expansion_ave=np.mean(expansion,1)
expansion_err=np.std(expansion,1)/np.sqrt(K)
    
plt.figure(figsize=(6, 6))

plt.errorbar(np.arange(2,num_T)*50,expansion_ave,yerr=expansion_err,linewidth=3)
plt.xlabel('Temperature (K)', fontsize=15)
plt.ylabel('Expansion coefficient (1/K)', fontsize=15)

plt.tight_layout()
plt.savefig('fig-chapter-12-expansion2.pdf')