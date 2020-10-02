clear; close all; font_size = 12;

% read in position data
load xf.txt;r=xf(752*64+1:end,1:3);

% some parameters to be set
N=64; % number of particles
L=[11,11,11]; % box lengths
L_times_pbc=L.*[1 1 1]; % boundary conditions
dt=2; % time step in fs

% some parameters automatically calculated
M=size(r,1)/N % number of frames

% unwrap the coordinates and get the velocities:
v=zeros((M-1)*N,3);
for m=1:M-1
    for n=1:N
        dr=r(m*N+n,:)-r((m-1)*N+n,:);
        dr=dr-round(dr./L).*L_times_pbc;
        r(m*N+n,:)=r((m-1)*N+n,:)+dr;
        v((m-1)*N+n,:)=dr/dt; % A/fs
    end
end

% get the VAC and DOS
v_all = zeros(N, 3, M-1);  % all the velocity data
for n = 1 : M-1           % put the data into a good form
   v_all(:, :, n) = v((n - 1) * N + 1 : n * N, :);
end
Nc=500;
omega=1:1:200;
nu=omega/2/pi;
[vac, dos] = find_vac_and_dos(v_all, Nc, dt/1000, omega);

% get the Diffusion coefficient
vac = vac*1.0e4; % from A^2/fs^2 to nm^2/ps^2
D = zeros(Nc, 1);
for n = 1 : Nc
   D(n) = sum(vac(1:n)) * dt/1000 / 3;
end
D = D - 0.5 * dt/1000 * vac(1);

% plot the results:
t = (0 : Nc - 1) * dt/1000;

figure;
plot(t + dt/1000, D, '-', 'linewidth', 2);
ylim([0,0.03]);
xlabel('correlation time (ps)', 'fontsize', font_size);
ylabel('D (nm^2/ps)', 'fontsize', font_size);
set(gca, 'fontsize', font_size)
grid on;

figure;
plot(t, vac, '-', 'linewidth', 2);
xlabel('Correlation Time (ps)', 'fontsize', font_size);
ylabel('VAC (nm^2/ps^2)', 'fontsize', font_size);
set(gca, 'fontsize', font_size);
axis tight;
grid on;

figure;
plot(nu, dos, '-', 'linewidth', 2);
xlabel('\nu (THz)', 'fontsize', font_size);
ylabel('PDOS (1/THz)', 'fontsize', font_size);
set(gca, 'fontsize', font_size);
axis tight;
grid on;

% should be approximately normalized to 1
normalization = sum(dos) * (nu(2) - nu(1))


