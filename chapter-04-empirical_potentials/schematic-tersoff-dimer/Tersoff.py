import matplotlib.pyplot as plt
import numpy as np

A = 1.83e3
B = 4.71e2
Lambda = 2.48
mu = 1.73
R = 2.7
S = 3.0

r = np.linspace(1.5, 4, 200)


def f_c(r):
    y = 0.5 * (1 + np.cos(np.pi * (r - R) / (S - R)))
    y[r < R] = 1
    y[r > S] = 0
    return y


def f_c_dev(r):
    y = 0.5 * np.pi / (S - R) * (-np.sin(np.pi * (r - R) / (S - R)))
    y[r < R] = 0
    y[r > S] = 0
    return y


def U_inner(r):
    return A * np.exp(-Lambda * r) - B * np.exp(-mu * r)


def U_inner_dev(r):
    return -Lambda * A * np.exp(-Lambda * r) + mu * B * np.exp(-mu * r)


def U(r):
    return f_c(r) * U_inner(r)


def F(r):
    f = f_c_dev(r) * U_inner(r) + f_c(r) * U_inner_dev(r)
    return -f


fig, ax1 = plt.subplots()
color = "tab:blue"
ax1.plot(r, U(r), color=color)
ax1.set_xlim(1.5, 3.5)
ax1.set_ylim(-3, 1)
ax1.set_xlabel("Distance (Å)")
ax1.set_ylabel("Energy (eV)", color=color)
plt.hlines(0, 1.5, 3.5, colors="k", zorder=-2, linestyles=":")
ax1.set_xticks([2, 2.5, 3])
ax1.set_yticks(range(-3, 2))

color = "tab:orange"
ax2 = ax1.twinx()
ax2.plot(r, F(r), color=color)
ax2.set_ylabel("Force (eV/Å)", color=color)
ax2.set_ylim(-15, 5)
ax2.set_yticks(range(-15, 6, 5))
ax2.vlines([R, S], [-15, -15], [5, 5], colors="grey", linestyles="--", zorder=-1)
ax2.text(R, -14, " R")
ax2.text(S, -14, " S")
plt.savefig("fig-chapter-8-tersoff.pdf")
