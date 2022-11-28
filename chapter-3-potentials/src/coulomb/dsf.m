function [energy,force]=dsf(num_atoms,box,r,q,rc,alpha)

energy=0;
force=zeros(3,num_atoms);
  
tic;
for i=1:num_atoms-1
    for j=i+1:num_atoms
        rij=r(:,j)-r(:,i);
        rij=rij-round(rij./box).*box;
        dij=norm(rij);
        if dij<rc
            Uij = erfc(alpha*dij )/dij - erfc(alpha * rc)/rc;
            Uij = Uij + erfc(alpha * rc) / (rc * rc) * (dij - rc);
            Uij = Uij * q(i)*q(j);
            energy=energy+Uij;
            f12=-(2/sqrt(pi)*alpha*dij*exp(-alpha*alpha*dij*dij)+erfc(alpha*dij))*rij/(dij*dij*dij);
            f12 = f12 + erfc(alpha * rc) / (rc * rc) * rij/dij;
            f12 = f12 * q(i)*q(j);
            force(:,i)=force(:,i)+f12;
            force(:,j)=force(:,j)-f12;
        end
    end
end
toc;

KC=14.399645;
energy=energy*KC;
force=force*KC;
end
