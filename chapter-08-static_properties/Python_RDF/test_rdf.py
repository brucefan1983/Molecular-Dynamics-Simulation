# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 21:01:27 2024

@author: nan xu (tamas@zju.edu.cn)
"""

import numpy as np
from matplotlib import pyplot as plt

def find_rdf(type1,type2,type,traj,box,pbc,Ng,rc):
    # determine some parameters
    N=traj.shape[0]    # number of particles
    dr=rc/Ng           # bin size
    N_type1=int(sum((type==type1))) # number of type1 atoms
    N_type2=int(sum((type==type2))) # number of type2 atoms
    volume = np.linalg.det(box)
    rho=N_type2/volume              #particle density for type 2

    # accumulate
    g=np.zeros((Ng,1))
    for n1 in range(1,N+1):         # sum over the atoms
        if type[n1-1] != type1:        # type1 is the center atom
            continue
        for n2 in range(1,N+1):     # loop over the atoms again
            if type[n2-1] != type2 or n1==n2: # type2 is the other atom 
                continue
            r12=traj[n2-1,:]-traj[n1-1,:]     # position difference vector
            r12=r12.T                         # column vector    
            r12=np.matmul(np.linalg.inv(box),r12) # transform to cubic box
            r12=r12-np.array(pbc)*np.round(r12) # mininum image convention
            r12=np.matmul(box,r12)              #transform back
            d12=np.sqrt((r12*r12).sum());       # distance
            if d12<rc:                           # there is a cutoff
                index=int(np.ceil(d12/dr))       # bin index
                g[index-1]=g[index-1]+1         # accumulate

    
    #normalize
    for n in range(1,Ng+1):
        g[n-1]=g[n-1]/N_type1              # average over the center atoms
        dV=4*np.pi*(dr*n)**2*dr            # volume of a spherical shell
        g[n-1]=g[n-1]/dV                   # now g is the local density
        g[n-1]=g[n-1]/rho                  # now g is the RDF
    return g

if __name__ == "__main__":
    # read in position data
    traj = np.loadtxt("r.xyz")
    N=132               # number of particles
    #your atom types
    # C H N O Ru Zr
    # 0 1 2 3 4  5
    type=traj[:N,0]     # atom types
    traj=traj[:,1:4]    # coordinates in units of Angstrom
    
    # choose the atom types you want to consider:
    type1=0             # C
    type2=4             # Ru
    
    # parameters from MD simulation
    pbc=[1,1,1]         # boundary conditions, 1 is periodic
    # a=5.6*4;
    # b=5.6*4;
    # c=5.6*4;
    # alpha=pi/2;
    # beta=pi/2;
    # gamma=pi/2;
    a=14.6599
    b=14.6518
    c=14.6567
    alpha=59.9849/180*np.pi
    beta=59.9923/180*np.pi
    gamma=59.9893/180*np.pi
    
    # choose the number of bins (number of data points in the figure)
    Ng=100         
    
    # determine other parameters automatically
    ax=a
    bx=b*np.cos(gamma)
    by=b*np.sin(gamma)
    cx=c*np.cos(beta)
    cy=(b*c*np.cos(alpha)-bx*cx)/by
    cz=np.sqrt(c*c-cx*cx-cy*cy)
    box=np.array([[ax,bx,cx],[0,by,cy],[0,0,cz]]) # box matrix
    volume=np.linalg.det(box)
    h1=volume/(a*b*np.sin(gamma))
    h2=volume/(b*c*np.sin(alpha))
    h3=volume/(c*a*np.sin(beta))
    rc=min([h1,h2,h3])/2  # the maximum cutoff radius that can be considered
    dr=rc/Ng              # bin size
    Ns=int(traj.shape[0]/N)    # number of frames
    
    # do the calculatigons
    g=np.zeros((Ng,1))    # The RDF to be calculated
    for n in range(1,Ns+1):
        r1=traj[((n-1)*N):(n*N),:]   # positions in one frame
        g=g+find_rdf(type1,type2,type,r1,box,pbc,Ng,rc) #sum over frames
    
    g=g/Ns  # time average in MD
    
    # plot the data
    r=np.linspace(1,Ng,Ng)*dr
    fig, axs = plt.subplots(1, 1)
    axs.plot(r,g,linewidth=3,marker="o",markerfacecolor='none', \
             markeredgewidth=2,markersize=7)
    axs.set_xlabel("r (Angstrom)",fontsize=15)
    axs.set_ylabel("g(r)",fontsize=15)
    axs.set_title("C-Ru",fontsize=15)
    axs.tick_params(axis='x', labelsize=15)
    axs.tick_params(axis='y', labelsize=15)
    axs.set_aspect(3/3.2, adjustable='box')
    plt.savefig("rdf.pdf", bbox_inches="tight")
    plt.show()