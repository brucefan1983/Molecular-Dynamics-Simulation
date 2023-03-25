clear; close all;
temperature=0.1:0.1:1; beta=1./temperature;
energy_theory=0.5+exp(-beta)./(1-exp(-beta));
C_theory=beta.^2.*exp(-beta)./(1-exp(-beta)).^2;

n_beads=32;
energy_32beads=zeros(1,length(temperature));
for n=1:length(temperature)
    disp(n)
    tic;
    energy_32beads(n)=get_energy(n_beads,beta(n));
    toc;
end

save('energy_32beads','energy_32beads');

load energy_4beads;
load energy_8beads;
load energy_16beads;
figure;
plot(1./beta,energy_4beads,'d-','linewidth',2); hold on;
plot(1./beta,energy_8beads,'s-','linewidth',2); hold on;
plot(1./beta,energy_16beads,'o-','linewidth',2); hold on;
plot(1./beta,energy_32beads,'x-','linewidth',2); hold on;
plot(1./beta,energy_theory,'^-','linewidth',2)
xlabel('Temperature');
ylabel('Energy');
legend('4 beads','8 beads','16 beads','32 beads','theory');
set(gca,'fontsize',18);


