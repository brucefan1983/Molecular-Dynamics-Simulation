import numpy as np

r0 = np.array([[0.0, 0.0, 0.5, 0.5],
               [0.0, 0.5, 0.0, 0.5],
               [0.0, 0.5, 0.5, 0.0]]).T

n0 = r0.shape[0]
nxyz = 6 * np.array([1, 1, 1])
N = nxyz[0] * nxyz[1] * nxyz[2] * n0
a = 5.385 * np.array([1, 1, 1])
box_length = a * nxyz

r = np.zeros((N, 3))
n = 0

for nx in range(nxyz[0]):
    for ny in range(nxyz[1]):
        for nz in range(nxyz[2]):
            for m in range(n0):
                n += 1
                r[n-1, :] = a * (np.array([nx, ny, nz]) + r0[m, :])

with open('xyz.in', 'w') as fid:
    fid.write(f'{N}\n')
    fid.write(f'{box_length[0]} {box_length[1]} {box_length[2]}\n')
    for n in range(N):
        fid.write(f'Ar {r[n, 0]} {r[n, 1]} {r[n, 2]} 40\n')
