import numpy as np
import matplotlib.pyplot as plt

data=[]
data.append(np.loadtxt('7A/thermo.out'))
data.append(np.loadtxt('9A/thermo.out'))
data.append(np.loadtxt('12A/thermo.out'))
data.append(np.loadtxt('15A/thermo.out'))

fluctuation_value = np.zeros(4)
for n in range(4):
    energy = np.sum(data[n][1:,1:3], axis = 1)
    fluctuation_value[n] = np.std((energy - np.mean(energy))/np.abs(energy))
    print(fluctuation_value[n])

plt.figure(figsize=(8, 7))

plt.semilogy(np.array((7,9,12,15)), fluctuation_value, '-o', linewidth=2,markersize=10)
plt.xlabel('Cutoff radius ($\mathrm{\AA}$)', fontsize=20)
plt.ylabel('Magnitude of relative energy fluctuations', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.tight_layout()
plt.savefig('fig-c02-cutoff.pdf')
