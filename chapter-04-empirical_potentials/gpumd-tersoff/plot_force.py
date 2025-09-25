import matplotlib.pyplot as plt
import numpy as np

force = np.loadtxt('force.out')
force_cpp = np.loadtxt('../cpp-tersoff-validation/force.out')

plt.figure()
plt.plot(force[0:64,:]-force_cpp[0:64,:], 'o')
plt.xlabel('atom index')
plt.ylabel('Force difference (eV/Å)')
plt.savefig('fig-chapter-8-tersoff-gpu-vs-cpu.pdf')
