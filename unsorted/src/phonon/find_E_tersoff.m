function energy=find_E_tersoff(n1,n2,L,pbc,r,type,NN,NL,para)
% may use type later
a=para(1);
b=para(2);
lambda=para(3);
mu=para(4);
beta=para(5);
n=para(6); 
c=para(7);
d=para(8);
h=para(9);
r1=para(10);
r2=para(11);
c2=c*c;
d2=d*d;
one_plus_c2overd2=1.0+c2/d2;
pi_factor=pi/(r2-r1);
minus_half_over_n=-0.5/n;
L_time_pbc=L.*pbc;
list=unique([n1,n2,NL(n1,1:NN(n1)),NL(n2,1:NN(n2))]);
energy = 0;
for i1=1:length(list)
    n1=list(i1);
    for i2=1:NN(n1)
        n2=NL(n1,i2);
        r12=r(n2,:)-r(n1,:);
        r12=r12-round(r12./L).*L_time_pbc; % minimum image convention
        d12=norm(r12);
        fr=a*exp(-lambda*d12);
        fa=b*exp(-mu*d12);
        if d12<r1
            fc=1.0;
        else
            fc=cos(pi_factor*(d12-r1))*0.5+0.5;
        end
        zeta=0.0;
        for i3=1:NN(n1)
            n3=NL(n1,i3);
            if n3==n2
                continue;
            end
            r13=r(n3,:)-r(n1,:);
            r13=r13-round(r13./L).*L_time_pbc; % minimum image convention
            d13=norm(r13);
            cos123=sum(r12.*r13)/(d12*d13);
            if d13<r1
                fc13=1.0;
            else
                fc13=cos(pi_factor*(d13-r1))*0.5+0.5;
            end
            g=one_plus_c2overd2-c2/(d2+(cos123-h)*(cos123-h));
            zeta=zeta+fc13*g;
        end
        b12=(1+(beta*zeta)^n)^minus_half_over_n;
        energy=energy+0.5*fc*(fr-b12*fa);
    end
end
