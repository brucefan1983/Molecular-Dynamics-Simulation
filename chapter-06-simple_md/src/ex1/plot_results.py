import numpy as np
import matplotlib.pyplot as plt

thermo = np.loadtxt('thermo.out')
timeStep = 5 / 1000  # ps
sampleInterval = 100
timeInterval = timeStep * sampleInterval
numData = thermo.shape[0]
time = np.arange(1, numData + 1) * timeInterval
totalEnergy = thermo[:, 1] + thermo[:, 2]
relativeEnergy = totalEnergy / np.mean(np.abs(totalEnergy))

plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(time, thermo[:, 1], '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Kinetic Energy (eV)', fontsize=15)
plt.title('(a)', fontsize=15)

plt.subplot(2, 2, 2)
plt.plot(time, thermo[:, 2], '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Potential Energy (eV)', fontsize=15)
plt.title('(b)', fontsize=15)

plt.subplot(2, 2, 3)
plt.plot(time, totalEnergy, '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Total Energy (eV)', fontsize=15)
plt.title('(c)', fontsize=15)

plt.subplot(2, 2, 4)
plt.plot(time[1:], relativeEnergy[1:], '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Relative Energy', fontsize=15)
plt.title('(d)', fontsize=15)

plt.tight_layout()
plt.savefig('fig-chapter-6-energy.pdf')

plt.figure()
plt.plot(time, thermo[:, 0], '-', linewidth=2)
plt.xlabel('Time (ps)', fontsize=15)
plt.ylabel('Temperature (K)', fontsize=15)

plt.tight_layout()
plt.savefig('fig-chapter-6-temperature.pdf')