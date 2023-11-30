clear; close all;
rng(1);
r0 = [0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0; ...
      0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0; ...
      0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5];
q0 = [1 1 1 1 -1 -1 -1 -1];
n0 = size(r0, 2);
nxyz = 4 * [1; 1; 1];
num_atoms = nxyz(1) * nxyz(2) * nxyz(3) * n0;
a = 6.5704 * [1; 1; 1];
box = a .* nxyz;

r = zeros(3,num_atoms);
q = zeros(num_atoms, 1);
n = 0;
for nx = 0 : nxyz(1) - 1
    for ny = 0 : nxyz(2) - 1
        for nz = 0 : nxyz(3) - 1
            for m = 1 : n0
                n = n + 1;
                r(:,n) = a .* ([nx;ny;nz] + r0(:,m))+rand;   
                q(n) = q0(m);
            end
        end
    end
end

figure;
plot3(r(1,:),r(2,:),r(3,:),'o');
axis equal;

kmax=3:9;
force=zeros(3,num_atoms,length(kmax));
for n=1:length(kmax)
    [energy,force(:,:,n)]=ewald(num_atoms,box,r,q,10,3/10,kmax(n));
    if n>1
        figure;
        plot(force(:,:,n).',force(:,:,n-1).','o');
    end
end


% finite difference checked
% force_fd=zeros(3,num_atoms);
% delta=1e-5;
% for n=1:1
%     for d=1:3
%         rp=r; rp(d,n)=rp(d,n)+delta;
%         rm=r; rm(d,n)=rm(d,n)-delta;
%         [energyp]=ewald(num_atoms,box,rp,q,rc,alpha,kmax);
%         [energym]=ewald(num_atoms,box,rm,q,rc,alpha,kmax);
%         force_fd(d,n)=(energym-energyp)/(2*delta);
%     end
% end
% force(:,1)-force_fd(:,1)

