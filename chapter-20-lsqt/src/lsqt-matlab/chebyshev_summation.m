function s=chebyshev_summation(M,C,E,E_max)
for m=0:M-1
    a=1/(M+1);
    g=(1-m*a)*cos(pi*m*a)+a*sin(pi*m*a)/tan(pi*a);
    C(m+1)=C(m+1)*g;
end
s=zeros(1,length(E));
for ne=1:length(E)
    T1=1;
    T2=E(ne);
    s(ne)=C(2)*T2;
    for m=3:M
        T_temp=T2;
        T2=2*E(ne)*T2-T1;
        T1=T_temp;
        s(ne)=s(ne)+C(m)*T2;
    end
    s(ne)=(s(ne)*2+C(1));
    s(ne)=s(ne)*2/(pi*sqrt(1-E(ne)*E(ne))*E_max);
end

