import matplotlib.pyplot as plt
import numpy as np

nu=np.loadtxt('nemd/nu.txt')
G=np.loadtxt('nemd/G.txt')
k=np.loadtxt('hnemd/k.txt')
nu=nu[:785]
G=G[:785]
k=k[:785]
lambda_omega=k/G
L=np.linspace(10,100000,10000);
kappa_L=np.zeros(10000)
for n in range(10000):
    kappa_L[n]=np.trapz(k/(1+lambda_omega/L[n]),nu)

plt.figure(figsize=(8, 8))
plt.subplot(2,1,1)
plt.semilogy(nu,lambda_omega/1000, '-', linewidth=2)
plt.title('(a)')
plt.ylim((0.01,10))
plt.xlabel('$\omega/2\pi$ (THz)', fontsize=15)
plt.ylabel('$\lambda$ ($\mu$m)', fontsize=15)
plt.tight_layout()

plt.subplot(2,1,2)
plt.semilogx(L/1000,kappa_L, '-', linewidth=2)
plt.title('(b)')
plt.xlabel('L ($\mu$m)', fontsize=15)
plt.ylabel('$\kappa$(L)', fontsize=15)
plt.tight_layout()

plt.savefig('fig-c9-kappa-lambda-omega.pdf')



