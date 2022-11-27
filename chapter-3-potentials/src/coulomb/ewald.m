function [energy,force]=ewald(num_atoms,box,r,q,rc,alpha,kmax)

V=box(1)*box(2)*box(3);
energy=0;
force=zeros(3,num_atoms);

tic;
for i=1:num_atoms-1
    for j=i+1:num_atoms
        rij=r(:,j)-r(:,i);
        rij=rij-round(rij./box).*box;
        dij=norm(rij);
        if dij<rc
            Uij=q(i)*q(j)*erfc(alpha*dij )/dij;
            energy=energy+Uij;
            f12=-q(i)*q(j)*(2/sqrt(pi)*alpha*dij*exp(-alpha*alpha*dij*dij)+erfc(alpha*dij))*rij/(dij*dij*dij);
            force(:,i)=force(:,i)+f12;
            force(:,j)=force(:,j)-f12;
        end
    end
end
toc;

tic;
num_kpoints=0;
k=zeros(3,10000);
for kx=-kmax:kmax
    for ky=-kmax:kmax
        for kz=-kmax:kmax
            k2=kx*kx+ky*ky+kz*kz;
            if k2<=kmax*kmax && k2>0
                num_kpoints=num_kpoints+1;
                k(:,num_kpoints)=2*pi*[kx;ky;kz]./box;
            end
        end
    end
end
k=k(:,1:num_kpoints);

G=zeros(1,num_kpoints);
for nk=1:num_kpoints
    ksq=sum(k(:,nk).*k(:,nk));
    G(nk)=(2*pi/V/ksq)*exp(-ksq/(4*alpha*alpha));
end

rho_per_atom=zeros(num_atoms,num_kpoints);
for na=1:num_atoms
    for nk=1:num_kpoints
        rho_per_atom(na,nk)=q(na)*exp(-1i*sum(k(:,nk).*r(:,na)));
    end
end
rho=sum(rho_per_atom);

energy=energy+real(sum(G.*(rho.*conj(rho))));
energy=energy-alpha*sum(q.*q)/sqrt(pi);
for na=1:num_atoms
    for nk=1:num_kpoints
        force(:,na)=force(:,na)+2*q(na)*G(nk)*k(:,nk)*imag(rho(nk)*exp(1i*sum(k(:,nk).*r(:,na))));
    end
end
toc;

KC=14.399645;
energy=energy*KC;
force=force*KC;
end
