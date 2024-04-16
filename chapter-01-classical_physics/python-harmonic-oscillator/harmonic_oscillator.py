import numpy as np
import matplotlib.pyplot as plt

m=1; k=1; dt=0.01; n_step=1000
v=0; x=1
v_vector = np.zeros(n_step)
x_vector = np.zeros(n_step)
for step in range(n_step):
    v = v + (dt/2) * (-k * x / m)
    x = x + dt * v
    v = v + (dt/2) * (-k * x / m)
    v_vector[step] = v
    x_vector[step] = x

plt.figure()
plt.plot(np.arange(1, n_step+1) * dt, x_vector, linewidth=2)
plt.plot(np.arange(1, n_step+1) * dt, v_vector, '--', linewidth=2)
plt.xlabel('time')
plt.ylabel('position or velocity')
plt.legend(['position', 'velocity'])
plt.savefig('fig-c01-position_and_momentum.pdf')

potential = 0.5 * k * x_vector**2
kinetic = 0.5 * m * v_vector**2

plt.figure()
plt.plot(np.arange(1, n_step+1) * dt, potential, linewidth=2)
plt.plot(np.arange(1, n_step+1) * dt, kinetic, '--', linewidth=2)
plt.plot(np.arange(1, n_step+1) * dt, potential + kinetic, '-.', linewidth=1)
plt.xlabel('time')
plt.ylabel('Energy')
plt.legend(['potential', 'kinetic', 'total'])
plt.savefig('fig-c01-energy_conservation.pdf')

plt.figure()

plt.figure()
plt.plot(x_vector, v_vector, '.', markersize=10)
plt.xlabel('position')
plt.ylabel('momentum')
plt.axis('equal')
plt.savefig('fig-c01-phase_space.pdf', dpi=200)
