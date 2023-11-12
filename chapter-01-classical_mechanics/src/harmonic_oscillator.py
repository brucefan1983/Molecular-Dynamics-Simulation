import numpy as np
import matplotlib.pyplot as plt

m=1; k=1; dt=0.01; n_step=1000
v=0; x=1
v_vector = np.zeros(n_step)
x_vector = np.zeros(n_step)
for step in range(1, n_step+1):
    v = v - (dt/2) * k * x
    x = x + dt * v
    v = v - (dt/2) * k * x
    v_vector[step-1] = v
    x_vector[step-1] = x

plt.figure()
plt.plot(np.arange(1, n_step+1) * dt, x_vector, linewidth=2)
plt.plot(np.arange(1, n_step+1) * dt, v_vector, '--', linewidth=2)
plt.xlabel('time')
plt.ylabel('position or velocity')
plt.legend(['position', 'velocity'])
plt.savefig('fig-chapter-1-position_and_momentum.pdf')

potential = 0.5 * k * x_vector**2
kinetic = 0.5 * m * v_vector**2

plt.figure()
plt.plot(np.arange(1, n_step+1) * dt, potential, linewidth=2)
plt.plot(np.arange(1, n_step+1) * dt, kinetic, '--', linewidth=2)
plt.plot(np.arange(1, n_step+1) * dt, potential + kinetic, '-.', linewidth=3)
plt.xlabel('time')
plt.ylabel('Energy')
plt.legend(['potential', 'kinetic', 'total'])
plt.savefig('fig-chapter-1-energy_conservation.pdf')

plt.figure()
plt.plot(x_vector, v_vector, '.', markersize=20)
plt.xlabel('position')
plt.ylabel('momentum')
plt.savefig('fig-chapter-1-phase_space.pdf', dpi=200)
