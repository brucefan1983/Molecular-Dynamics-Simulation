# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 19:58:00 2024

@author: nan xu (tamas@zju.edu.cn)
"""

import numpy as np
import time
from matplotlib import pyplot as plt

def get_energy(n_beads,beta):
    #core function
    hbar=1; m=1; _lambda=1; dt=0.1; tau_T=100
    omega_n=n_beads/beta/hbar; n_step=1000000; n_step_pimd=1000000;
    cayley=True  # cayley is much more stable
    C=np.zeros((n_beads,n_beads))
    for j in range(0,n_beads):
        for k in range(0,n_beads):
            if k==0:
                C[j,k]=np.sqrt(1/n_beads)
            elif k<=(n_beads/2-1):
                C[j,k]=np.sqrt(2/n_beads)*np.cos(2*np.pi*j*k/n_beads)
            elif k==n_beads/2:
                C[j,k]=np.sqrt(1/n_beads)*(-1)**j
            else:
                C[j,k]=np.sqrt(2/n_beads)*np.sin(2*np.pi*j*k/n_beads)
    p=np.zeros(n_beads); q=np.ones(n_beads) 
    energy=np.zeros((n_step,1))
    for step in range(1,n_step+1):   
        p_normal=np.matmul(p,C)
        c1=np.exp(-dt*omega_n*np.sin(np.linspace(0,n_beads-1,n_beads)*np.pi/n_beads))
        if step<=n_step_pimd:
            c1[0]=np.exp(-1/2/tau_T)
        c2=np.sqrt(1-c1**2)
        p_normal=c1*p_normal+np.sqrt(n_beads*m/beta)*c2* np.random.standard_normal(n_beads)  
        p=(np.matmul(C,p_normal.T)).T
        p=p-(dt/2)*m*_lambda*_lambda*q
        p_normal=np.matmul(p,C)
        q_normal=np.matmul(q,C) 
        for k in range(0,n_beads):
            omega_k=2*omega_n*np.sin(k*np.pi/n_beads) 
            if k==0:
                q_normal[k]=(dt/m)*p_normal[k]+q_normal[k]
            else:
                c=np.cos(omega_k*dt); s=np.sin(omega_k*dt)
                if cayley:
                    c=(1-(omega_k*dt/2)**2)/(1+(omega_k*dt/2)**2)
                    s=omega_k*dt/(1+(omega_k*dt/2)**2)
                p_temp=p_normal[k]
                q_temp=q_normal[k]
                p_normal[k]=c*p_temp-m*omega_k*s*q_temp
                q_normal[k]=(1/m/omega_k)*s*p_temp+c*q_temp
        p=(np.matmul(C,p_normal.T)).T; q=(np.matmul(C,q_normal.T)).T
        p=p-(dt/2)*m*_lambda*_lambda*q
        p_normal=np.matmul(p,C)    
        c1=np.exp(-dt*omega_n*np.sin(np.linspace(0,n_beads-1,n_beads)*np.pi/n_beads))
        if step<=n_step_pimd:
            c1[0]=np.exp(-1/2/tau_T)
        c2=np.sqrt(1-c1**2)
        p_normal=c1*p_normal+np.sqrt(n_beads*m/beta)*c2*np.random.standard_normal(n_beads,) 
        p=(np.matmul(C,p_normal.T)).T
        q_ave=np.mean(q)
        kinetic_energy=0.5/beta+0.5*m*_lambda*_lambda*np.mean((q-q_ave)*q)
        potential_energy=0.5*m*_lambda*_lambda*np.mean(q**2)
        energy[step-1]=kinetic_energy+potential_energy
    energy=np.mean(energy[int(len(energy)/2):])
    return energy

if __name__ == "__main__":
    temperature=np.linspace(0.1,1,10)
    beta=1./temperature
    energy_theory=0.5+np.exp(-beta)/(1-np.exp(-beta))
    C_theory=beta**2*np.exp(-beta)/(1-np.exp(-beta))**2
    results = []
    for n_beads in [4,8,16,32]:
        energy_n_beads=np.zeros(len(temperature))
        for n in range(1,len(temperature)+1):
            print(n)
            start_time = time.time()
            energy_n_beads[n-1]=get_energy(n_beads,beta[n-1])
            print("Consumed time: %f seconds." %(time.time()-start_time))
        results.append(energy_n_beads)
        np.save('energy_%dbeads' %(n_beads),energy_n_beads)

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
