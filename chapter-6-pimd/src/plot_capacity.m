clear; close all;
temperature=0.1:0.1:1; beta=1./temperature;
energy_quantum=0.5+exp(-beta)./(1-exp(-beta));
energy_classical=temperature;
load energy_4beads;
load energy_8beads;
load energy_16beads;
load energy_32beads;
figure;
plot(1./beta,energy_4beads,'d-','linewidth',2); hold on;
plot(1./beta,energy_8beads,'s-','linewidth',2); hold on;
plot(1./beta,energy_16beads,'o-','linewidth',2); hold on;
plot(1./beta,energy_32beads,'^-','linewidth',2); hold on;
plot(1./beta,energy_classical,'--','linewidth',2)
plot(1./beta,energy_quantum,'-','linewidth',2)
xlabel('Temperature');
ylabel('Energy');
legend('4 beads','8 beads','16 beads','32 beads','Classical Theory','Quantum Theory');
set(gca,'fontsize',18);


