function MSD=find_msd(M,E_max,dt,E,H,V,phi)
phix=phi*0;
MSD=zeros(length(dt),length(E));
for nt=1:length(dt)
    phix=evolve(H,dt(nt),1,phix);
    phix=phix+evolvex(H,V,dt(nt),phi);
    phi=evolve(H,dt(nt),1,phi);
    Cn=find_moments(M,H,phix,phix);
    MSD(nt,:)=chebyshev_summation(M,Cn,E,E_max);
end
