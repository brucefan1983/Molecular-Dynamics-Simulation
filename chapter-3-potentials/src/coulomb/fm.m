function [energy,force]=fm(num_atoms,box,r,q,rc)

energy=0;
force=zeros(3,num_atoms);

a=[-0.165477570871E-03 -0.547587014180E-04
    0.288823451703E-03 -0.777061023641E-05
    -0.122247561247E-03 0.649658451885E-04
    0.963712701767E-05 -0.417496679930E-04
    0.251954672874E-06 0.926924324623E-05
    -0.735796273353E-07 -0.107542095070E-05
    0.353601771929E-08 0.710779773104E-07
    -0.525765995765E-10 -0.261040455982E-08
    0.0 0.439931939572E-10
    0.0 -0.422656965444E-14
    0.0 -0.656184357691E-14];

tic;
for i=1:num_atoms-1
    for j=i+1:num_atoms
        rij=r(:,j)-r(:,i);
        rij=rij-round(rij./box).*box;
        dij=norm(rij);
        if dij<rc
            nij=rij/dij;
            dij=dij/0.529177210903; % atomic
            f12=1/dij/dij;
            if rc<11
                for k=1:8
                    f12=f12+a(k,1)*dij^(k-1);
                end
            else
                for k=1:11
                    f12=f12+a(k,2)*dij^(k-1);
                end
            end
            
            f12=-f12*q(i)*q(j)*nij;
            f12=f12*27.211386245988/0.529177210903; %eV/Bohr to eV/A
            force(:,i)=force(:,i)+f12;
            force(:,j)=force(:,j)-f12;
        end
    end
end
toc;

%KC=14.399645;
%energy=energy*KC;
%force=force*KC;
end
