function energy=find_E(n1,n2,L,pbc,r,type,NN,NL,pot)
if pot.model==0
    energy=find_E_tersoff(n1,n2,L,pbc,r,type,NN,NL,pot.para);
else
    energy=0;
end
