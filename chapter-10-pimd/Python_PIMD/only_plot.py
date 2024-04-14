# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 19:58:00 2024

@author: nan xu (tamas@zju.edu.cn)
"""

import numpy as np
from matplotlib import pyplot as plt


if __name__ == "__main__":
    temperature=np.linspace(0.1,1,10)
    beta=1./temperature
    energy_theory=0.5+np.exp(-beta)/(1-np.exp(-beta))

    energy_4beads  = np.load("energy_4beads.npy")
    energy_8beads  = np.load("energy_8beads.npy")
    energy_16beads = np.load("energy_16beads.npy")
    energy_32beads = np.load("energy_32beads.npy")

    fig, axs = plt.subplots(1, 1)
    axs.plot(1./beta,energy_4beads,linewidth=3,marker="d",markerfacecolor='none', \
             markeredgewidth=2,markersize=7,label='4 beads')
    axs.plot(1./beta,energy_8beads,linewidth=3,marker="s",markerfacecolor='none', \
             markeredgewidth=2,markersize=7,label='8 beads')
    axs.plot(1./beta,energy_16beads,linewidth=3,marker="o",markerfacecolor='none', \
             markeredgewidth=2,markersize=7,label='16 beads')     
    axs.plot(1./beta,energy_32beads,linewidth=3,marker="x",markerfacecolor='none', \
             markeredgewidth=2,markersize=7,label='32 beads')           
    axs.plot(1./beta,energy_theory,linewidth=3,marker="^",markerfacecolor='none', \
             markeredgewidth=2,markersize=7,label='theory')  
    axs.legend(loc="best")
    axs.set_xlabel("Temperature",fontsize=15)
    axs.set_ylabel("Energy",fontsize=15)
    axs.tick_params(axis='x', labelsize=15)
    axs.tick_params(axis='y', labelsize=15)
    axs.set_aspect(3/3.2, adjustable='box')
    plt.savefig("pimd.pdf", bbox_inches="tight")
    plt.show()
