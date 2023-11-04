function phi=create_state(N)
random_phase=rand(N,1)*2*pi;
phi=cos(random_phase)+sin(random_phase)*1i;   
phi=phi/norm(phi);