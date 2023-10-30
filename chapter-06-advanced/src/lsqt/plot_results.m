clear; close all; font_size=8;

% From GPUMD input
E_min=-8.1; % eV
E_max=8.1; % eV
Ne=10001; % number of energy points
Nt=1000; % number of time points
dt=0.1; % fs
N_atom=187200; % number of atoms
volume=885.424372829*766.8*3.35; % A^3

% calculate energy
E=linspace(E_min,E_max,Ne); % eV
E=E-(-0.12); % based on DOS, the Fermi energy is at -0.12 eV
N_065eV=floor(5001+(0.65-0.12)/1.62e-3);


% averaged data for LSQT output
load dos_ave.mat;
load sigma_ave.mat;
load velocity_ave.mat;

% calculate transport distribution function (TDF)
tdf=mean(sigma_ave(end/5+1:end,:)); %S/m

% calculate the transport properties at 300 K
T=300; % K
kT=8.617343e-5*T; % eV
mu=-1:0.01:1; % chemical potential points
num_mu=length(mu);
sigma_300K=zeros(num_mu,1);
Seebeck_300K=zeros(num_mu,1);
kappa_300K=zeros(num_mu,1);
ZT_300K=zeros(num_mu,1);
kappa_ph=6.7; % W/mK  (from HNEMD)
for n=1:num_mu
    mu_n=mu(n);
    x=(E-mu_n)/kT;
    dfdE_factor = (1/kT).*exp(x)./(exp(x)+1)./(exp(x)+1);
    X0=trapz(E,tdf.*dfdE_factor);
    X1=trapz(E,E.*tdf.*dfdE_factor);
    X2=trapz(E,E.*E.*tdf.*dfdE_factor);
    sigma_300K(n)=X0;
    Seebeck_300K(n)=-(1/T)*(X1./X0-mu_n); % V/K
    kappa_300K(n)=(X2-X1.^2./X0)/T; % W/mK
    ZT_300K=Seebeck_300K.^2.*sigma_300K*T./(kappa_300K+kappa_ph);
end

% plot the results
figure; 

subplot(3,2,1);
plot([0,(1:Nt)*dt],[0;sigma_ave(:,N_065eV)]*1e-6,'-','linewidth',2);hold on;
xlabel('t (fs)');
ylabel('\Sigma (10^6 S/m)');
set(gca,'fontsize',font_size,'linewidth',1,'xtick',0:20:100);
title('(a)')

subplot(3,2,2);
yyaxis left;
plot(E,tdf*1e-6,'-','linewidth',2); hold on;
ylabel('\Sigma (10^6 S/m)');
ylim([0,4])
yyaxis right;
plot(E,dos_ave,'linewidth',2);hold on;
ylabel('DOS (states/eV/atom)');
ylim([0,0.15])
xlim([-1,1])
xlabel('E (eV)');
set(gca,'fontsize',font_size,'xtick',-1:0.5:1,'linewidth',1);
title('(b)')

subplot(3,2,3)
plot(mu,sigma_300K*1e-6,'-','linewidth',2); hold on;
ylabel('\sigma (10^6 S/m)');
ylim([0,3])
xlim([-1,1])
xlabel('\mu (eV)');
set(gca,'fontsize',font_size,'xtick',-1:0.5:1,'linewidth',1);
title('(c)')

subplot(3,2,4)
plot(mu,Seebeck_300K*1e3,'-','linewidth',2); hold on;
ylabel('S (10^{-3}V/K)');
%ylim([0,4])
xlim([-1,1])
xlabel('\mu (eV)');
set(gca,'fontsize',font_size,'xtick',-1:0.5:1,'linewidth',1);
title('(d)')

subplot(3,2,5)
plot(mu,kappa_300K,'-','linewidth',2); hold on;
ylabel('\kappa_{el} (W/mK)');
ylim([0,15])
xlim([-1,1])
xlabel('\mu (eV)');
set(gca,'fontsize',font_size,'xtick',-1:0.5:1,'linewidth',1);
title('(e)')

subplot(3,2,6)
plot(mu,ZT_300K,'-','linewidth',2); hold on;
ylabel('zT');
ylim([0,0.5])
xlim([-1,1])
xlabel('\mu (eV)');
set(gca,'fontsize',font_size,'xtick',-1:0.5:1,'linewidth',1);
title('(f)')

