clear; close all;

% potential parameter
pot.model=0; % tersoff
pot.para=[1393.6,430,3.5333,2.2407,1.5724e-7,0.72751,38049,4.3484,-0.93,1.8,2.1];

% r-space data
a0=2.46; % lattice constant
axyz=[a0,a0*sqrt(3),3.35];
r0=a0*[1/2,0,0;0,sqrt(3)/6,0;0,sqrt(3)/2,0;1/2,sqrt(3)*2/3,0];
basis_label0=[2,1,2,1];
type0=[1,1,1,1];
nxyz=[3,3,1];
pbc=[1,1,0];
[r,basis_label,type,L]=find_r(r0,basis_label0,type0,nxyz,axyz);

% k_space data
a=[1/2,1/2,0;-sqrt(3)/2,sqrt(3)/2,0;0,0,1]*a0;
special_k=[0,0,0;1/2,0,0;1/3,1/3,0;0,0,0];
name_special_k={'$\Gamma$','M','K','$\Gamma$'};
Nk=100; % number of k points between two special ones
[K,k_norm]=find_k(Nk,special_k.',a);

% specify the primitive cell atoms
basis_atom=[2,3]; % A and B
basis_mass=[12,12];

% get the results
[nu]=matlady(pot,2.0,2.6,basis_atom,basis_mass,r,basis_label,type,L,pbc,K);
plot_dispersion(Nk,k_norm/(2*pi/a0),name_special_k,nu);
