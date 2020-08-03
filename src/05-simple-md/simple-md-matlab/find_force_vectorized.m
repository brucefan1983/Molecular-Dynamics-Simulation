function [f,U]=find_force_vectorized(N,D,NN,NL,L,pbc,r)
EPSILON=1.032e-2;% in units of eV (only for Argon)
SIGMA=3.405; % in units of Angstrom (only for Argon)
sigma_6=SIGMA^6;sigma_12=SIGMA^12;L_times_pbc=L.*pbc;
U=0; % initialize the total potential energy of the system
f=zeros(N,D); % initialize the total forces on each atom
for n1=1:N % loop over the atoms
    m=1:NN(n1);
    n2=NL(n1,m);% all the neighbors of n1
    r12=zeros(NN(n1),D);
    for d=1:D
        r12(:,d)=r(n2,d)-r(n1,d);
        r12(:,d)=r12(:,d)-round(r12(:,d)/L(d)).*L_times_pbc(d); %minimum image convention
    end
    d12_square=sum(r12.*r12,2);
    d_6=d12_square.^3;d_8=d12_square.*d_6;d_12=d_6.*d_6;d_14=d_6.*d_8;
    f12=(sigma_6./d_8-2.0*sigma_12./d_14)*(24.0*EPSILON);
    f(n1,:)=f(n1,:)+sum([f12,f12,f12].*r12);
    U=U+2.0*EPSILON*sum(sigma_12./d_12-sigma_6./d_6); % avoid double counting
end




