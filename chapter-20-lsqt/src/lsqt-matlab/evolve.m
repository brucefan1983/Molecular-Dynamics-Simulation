function [phi_o]=evolve(H,dt,sign,phi_i)
phi_0=phi_i;
phi_1=H*phi_i;
j0=besselj(0,dt);
j1=besselj(1,dt);
phi_o=j0*phi_0+2*(-1i*sign)*j1*phi_1;
for m=2:100000000
    jm=besselj(m,dt);
    if abs(jm)<1.0e-15;break;end
    phi_2=2*H*phi_1-phi_0;
    phi_o=phi_o+2*(-1i*sign)^m*jm*phi_2;
    phi_0=phi_1;  
    phi_1=phi_2;
end
