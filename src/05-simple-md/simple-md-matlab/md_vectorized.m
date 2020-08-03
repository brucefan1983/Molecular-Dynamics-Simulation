function [E]=md_vectorized(D,n0,nxyz,a,pbc,Ne,Np,Ns,rc,dt,T)
K_B=8.625e-5; % Boltzmann's constant in my unit system  
TIME_UNIT_CONVERSION=10.18; % from fs to my unit system
N=n0*nxyz(1)*nxyz(2)*nxyz(3); % number of atoms
L=a.*nxyz; % box size (Angstrom)
dt=dt/TIME_UNIT_CONVERSION; % time step in my unit system
m=ones(N,1)*40; % mass for Argon atom in my unit system
r=initialize_position(N,D,n0,nxyz,a); % intial positions
v=initialize_velocity(K_B,N,D,T,m); % initial velocities
[NN,NL]=find_neighbor(N,L,pbc,rc,r); % fixed neighbor list
[f]=find_force_vectorized(N,D,NN,NL,L,pbc,r); % initial forces
E=zeros(Np/Ns,3); % energy data to be computed
for step=1:(Ne+Np) % time-evolution started
    for d=1:D % step 1 of Velocity-Verlet
        v(:,d)=v(:,d)+(f(:,d)./m)*(dt*0.5); 
        r(:,d)=r(:,d)+v(:,d)*dt; 
    end
    [f,U]=find_force_vectorized(N,D,NN,NL,L,pbc,r); % update forces
    for d=1:D % step 2 of Velocity-Verlet
        v(:,d)=v(:,d)+(f(:,d)./m)*(dt*0.5);
    end 
    if step<=Ne; % control temperature in the equilibration stage
        v=v*sqrt(T*D*K_B*N/sum(m.*sum(v.^2,2))); % scale velocity
    elseif mod(step,Ns)==0 % measure in the production stage
        E((step-Ne)/Ns,1)=U/N; % potential (per atom)
        E((step-Ne)/Ns,2)=0.5*sum(m.*sum(v.^2,2))/N; % kinetic energy
    end
end % time-evolution completed
E(:,3)=E(:,1)+E(:,2); % total enegy (per atom)
