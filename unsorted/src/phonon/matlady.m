function [nu]=matlady(pot,r1,r2,basis_atom,basis_mass,r,basis_label,type,L,pbc,K)
[NN1,NL1]=find_neighbor(r,1,r1,L,pbc);       % for force
[NN2,NL2]=find_neighbor(r,-1,r2,L,pbc);      % for force constant
tic;[Phi,R]=find_phi_all(pot,L,pbc,r,type,basis_atom,NN1,NL1,NN2,NL2);toc;
tic;[nu]=find_nu(NN2,NL2,Phi,R,basis_label,basis_atom,basis_mass,K);toc;
