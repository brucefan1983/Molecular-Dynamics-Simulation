function [Phi,R]=find_phi_all(pot,L,pbc,r,type,basis_atom,NN1,NL1,NN2,NL2) 
Phi=cell(length(basis_atom),1); % force constants 
R=cell(length(basis_atom),1); % relative positions
L_time_pbc=L.*pbc;
for m=1:length(basis_atom)
    atom=basis_atom(m); % atom-m in the basis
    Phi{m}=zeros(3,3,NN2(atom));
    R{m}=zeros(NN2(atom),3);
    for n=1:NN2(atom) % only consider nonzero force constants
        Phi{m}(:,:,n)=find_phi_one(atom,NL2(atom,n),r,type,NN1,NL1,L,pbc,pot);
        R{m}(n,:)=r(NL2(atom,n),:)-r(atom,:);
        R{m}(n,:)=R{m}(n,:)-round(R{m}(n,:)./L).*L_time_pbc;
    end
end
