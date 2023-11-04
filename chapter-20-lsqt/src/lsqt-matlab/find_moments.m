function C=find_moments(M,H,phi_left,phi_right)
C=zeros(1,M);
phi_1=phi_right;
C(1)=phi_left'*phi_1;
phi_2=H*phi_1;
C(2)=phi_left'*phi_2;
for m=3:M
    temp=phi_2;
    phi_2=2.0*H*phi_2-phi_1;
    phi_1=temp;
    C(m)=phi_left'*phi_2;
end
C=real(C);
