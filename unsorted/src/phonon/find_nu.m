function [nu]=find_nu(NN2,NL2,Phi,R,type,basis_atom,basis_mass,K)
Nb=size(basis_atom,2);
NK=size(K,2);
nu=zeros(3*Nb,NK);
for nk=1:NK
    k=K(:,nk);
    DK=zeros(Nb*3,Nb*3);
    for m=1:Nb
        atom_m=basis_atom(m);
        type_m=type(atom_m);
        mass_m=basis_mass(type_m);
        index_m=(type_m-1)*3+1:type_m*3;
        for n=1:NN2(basis_atom(m))
            atom_n=NL2(atom_m,n);
            type_n=type(atom_n);
            mass_n=basis_mass(type_n);
            index_n=(type_n-1)*3+1:type_n*3;
            DK(index_m,index_n)=DK(index_m,index_n)+Phi{m}(:,:,n)...
                *exp(1i*dot(k,R{m}(n,:)))/sqrt(mass_m*mass_n);
        end
    end
    if norm(DK-DK')>0.00001
        disp(norm(DK-DK'))
    end
    nu(:,nk)=sqrt(eig(DK));
end
nu=nu*1000/10.18/2/pi; % in units of THz now
