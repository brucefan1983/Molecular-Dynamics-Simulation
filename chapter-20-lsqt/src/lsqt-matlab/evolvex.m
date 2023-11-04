function [phi_o]=evolvex(H,V,dt,phi_i)
phi_0=phi_i;
phix_0=0;
phi_1=H*phi_0;
phix_1=1i*(V*phi_i);
phi_o=2*(-1i)*besselj(1,dt)*phix_1;
for m=2:100000000
    jm=besselj(m,dt);
    if abs(jm)<1.0e-15;break;end
    phi_2=2.0*H*phi_1-phi_0;
    phix_2=2*1i*(V*phi_1)+2*H*phix_1-phix_0;
    phi_o=phi_o+2*(-1i)^m*jm*phix_2;
    phi_0=phi_1;  
    phi_1=phi_2;
    phix_0=phix_1; 
    phix_1=phix_2;
end
