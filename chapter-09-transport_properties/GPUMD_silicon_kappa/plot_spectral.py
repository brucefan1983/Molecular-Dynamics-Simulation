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

plt.figure(figsize=(12, 6))
plt.subplot(1,2,1)
plt.semilogy(nu,lambda_omega/1000, '-', linewidth=3)
plt.title('(a)', fontsize=20)
plt.ylim((0.01,10))
plt.xlabel('$\omega/2\pi$ (THz)', fontsize=25)
plt.ylabel('$\lambda$ ($\mu$m)', fontsize=25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()

plt.subplot(1,2,2)
plt.semilogx(L/1000,kappa_L, '-', linewidth=3)
plt.title('(b)', fontsize=20)
plt.xlabel('L ($\mu$m)', fontsize=25)
plt.ylabel('$\kappa$(L)', fontsize=25)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()

plt.savefig('fig-c09-kappa-lambda-omega.pdf')



