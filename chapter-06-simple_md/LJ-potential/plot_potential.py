import matplotlib.pyplot as plt
import numpy as np


def U(r, rc):
    ret = 4 * (1 / r**12 - 1 / r**6)
    ret[r > rc] = 0
    return ret


def F(r, rc):
    ret = -4 * (-12 / r**13 + 6 / r**7)
    ret[r > rc] = 0
    return ret


rc = 4
r = np.linspace(1, 4, 500)

fig, ax1 = plt.subplots()

color = "tab:blue"
ax1.plot(r, U(r, rc), color=color)
ax1.hlines(0, r.min(), r.max(), linestyles="dashed", colors="k", alpha=0.5, zorder=-1)
ax1.set_xlabel(r"Distance ($\sigma$)")
ax1.set_ylabel(r"Energy ($\epsilon$)", color=color)
ax1.set_ylim(-1.2, 1.2)

color = "tab:orange"
ax2 = ax1.twinx()
ax2.plot(r, F(r, rc), color=color)
ax2.set_ylabel(r"Force ($\epsilon$/$\sigma$)", color=color)
ax2.set_ylim(-25, 25)

plt.savefig("fig-chapter-6-lj.pdf")
