clear; close all;
rng(1);
r0 = [0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.0; ...
      0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0; ...
      0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5];
q0 = [1 1 1 1 -1 -1 -1 -1];
n0 = size(r0, 2);
nxyz = 5 * [1; 1; 1];
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

[energy_ewald,force_ewald]=ewald(num_atoms,box,r,q,10,3/10,8);
[energy_dsf10A,force_dsf10A]=dsf(num_atoms,box,r,q,10,0.12);
[energy_fm10A,force_fm10A]=fm(num_atoms,box,r,q,10);

[energy_dsf13A,force_dsf13A]=dsf(num_atoms,box,r,q,13,0.1);
[energy_fm12A,force_fm12A]=fm(num_atoms,box,r,q,12);

figure;
plot(force_ewald.',force_dsf10A.','ro'); hold on;
plot(force_ewald.',force_fm10A.','bx'); hold on;
plot(linspace(-3,3,100),linspace(-3,3,100),'-');
xlabel('Ewald');
ylabel('DSF or FM');
legend('DSF-10A','FM-10A');
set(gca,'fontsize',15);

figure;
plot(force_ewald.',force_dsf13A.','ro'); hold on;
plot(force_ewald.',force_fm12A.','bx'); hold on;
plot(linspace(-3,3,100),linspace(-3,3,100),'-');
xlabel('Ewald');
ylabel('DSF or FM');
legend('DSF-12A','FM-12A');
set(gca,'fontsize',15);

sqrt(mean(mean((force_ewald-force_dsf10A).^2)))
sqrt(mean(mean((force_ewald-force_fm10A).^2)))
sqrt(mean(mean((force_ewald-force_dsf13A).^2)))
sqrt(mean(mean((force_ewald-force_fm12A).^2)))


