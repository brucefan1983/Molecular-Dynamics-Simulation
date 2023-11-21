function [H]=find_H1(pot,n1,n2,r0,NN,NL,L)
H=zeros(3,3);
delta=0.001; % small number used in finite difference
factor=1/delta/delta/4;
for alpha=1:3 % 3 directions for n1
    for beta=1:3 % 3 directions for n2
        % copy the positions first
        rpp=r0;rmm=r0;rpm=r0;rmp=r0;
        % make small shifts
        rpp(n1,alpha)=rpp(n1,alpha)+delta; % plus
        rpp(n2,beta)=rpp(n2,beta)+delta; % plus
        rmm(n1,alpha)=rmm(n1,alpha)-delta; % minus
        rmm(n2,beta)=rmm(n2,beta)-delta; % minus
        rpm(n1,alpha)=rpm(n1,alpha)+delta; % plus
        rpm(n2,beta)=rpm(n2,beta)-delta; % minus
        rmp(n1,alpha)=rmp(n1,alpha)-delta; % minus
        rmp(n2,beta)=rmp(n2,beta)+delta; % plus
        % calculate the energies
        epp=find_E(n1,n2,L,rpp,NN,NL,pot.para);
        emm=find_E(n1,n2,L,rmm,NN,NL,pot.para);
        epm=find_E(n1,n2,L,rpm,NN,NL,pot.para);
        emp=find_E(n1,n2,L,rmp,NN,NL,pot.para);
        % use finite difference to obtain the force constants
        H(alpha,beta)=(epp+emm-epm-emp)*factor;
    end
end