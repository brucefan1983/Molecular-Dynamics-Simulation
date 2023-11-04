function [Phi]=find_phi_one(n1,n2,r0,type,NN,NL,L,pbc,pot)
Phi=zeros(3,3);
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
        epp=find_E(n1,n2,L,pbc,rpp,type,NN,NL,pot);
        emm=find_E(n1,n2,L,pbc,rmm,type,NN,NL,pot);
        epm=find_E(n1,n2,L,pbc,rpm,type,NN,NL,pot);
        emp=find_E(n1,n2,L,pbc,rmp,type,NN,NL,pot);
        % use finite difference to obtain the force constants
        Phi(alpha,beta)=(epp+emm-epm-emp)*factor;
    end
end
