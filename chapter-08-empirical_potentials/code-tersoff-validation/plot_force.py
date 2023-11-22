import matplotlib.pyplot as plt
import numpy as np

force = np.loadtxt('force.out')


plt.figure()
plt.plot(force[64:,:]-force[0:64,:], 'o')
#plt.plot(force[64:,:], 'x')
plt.xlabel('atom index')
plt.ylabel('Force difference (eV/Å)')
plt.savefig('fig-chapter-8-tersoff-validation.pdf')
