function VAC=find_vac(M,E_max,dt,E,H,V,phi)
phi_left=phi;
phi_right=V*phi;
VAC=zeros(length(dt),length(E));
for nt=1:length(dt)
    C=find_moments(M,H,V*phi_left,phi_right);
    VAC(nt,:)=chebyshev_summation(M,C,E,E_max);
    phi_left=evolve(H,dt(nt),-1,phi_left); 
    phi_right=evolve(H,dt(nt),-1,phi_right);
end
