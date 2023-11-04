clear; close all;

% potential parameter
pot.model=0;
pot.para=[1.8308e3,471.18,2.4799,1.7322,1.1e-6,0.78734,1.0039e5,...
    16.217,-0.59825,2.7,3];

% r-space data
a0=5.432; % lattice constant
axyz=a0*[1,1,1];
r0=[0.0, 0.0, 0.5, 0.5, 0.25, 0.25, 0.75, 0.75; ...
    0.0, 0.5, 0.0, 0.5, 0.25, 0.75, 0.25, 0.75; ...
    0.0, 0.5, 0.5, 0.0, 0.25, 0.75, 0.75, 0.25].';
r0=r0*a0;
basis_label0=[1,1,1,1,2,2,2,2];
type0=[1,1,1,1,1,1,1,1];
nxyz=[3,3,3];
pbc=[1,1,1];
[r,basis_label,type,L]=find_r(r0,basis_label0,type0,nxyz,axyz);

% k-space data
a=[0 1 1;1 0 1;1 1 0]*a0/2;
special_k=[0,0,0;    1/2,0,1/2;    % Gamma -> X
    1/2,0,1/2;       5/8,1/4,5/8;  % X -> U=K
    3/8,3/8,3/4;     0,0,0;        % K -> Gamma
    0,0,0;           1/2,1/2,1/2]; % Gamma -> L
name_special_k={'$\Gamma$','X','U=K','$\Gamma$','L'};
Nk=100; % number of k points between two special ones
[K,k_norm]=find_k(Nk,special_k.',a);

% specify the primitive cell atoms
basis_atom=[1,5]; % A and B
basis_mass=[28,28];

% get the results
[nu]=matlady(pot,3.0,4.0,basis_atom,basis_mass,r,basis_label,type,L,pbc,K);
plot_dispersion(Nk,k_norm/(2*pi/a0),name_special_k,nu);
