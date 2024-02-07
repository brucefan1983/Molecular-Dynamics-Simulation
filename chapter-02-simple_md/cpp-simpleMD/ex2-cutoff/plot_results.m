clear; close all;
load 7A/thermo.out; energy_7A = sum(thermo(2:end,2:3),2);
load 9A/thermo.out; energy_9A = sum(thermo(2:end,2:3),2);
load 12A/thermo.out; energy_12A = sum(thermo(2:end,2:3),2);
load 15A/thermo.out; energy_15A = sum(thermo(2:end,2:3),2);

fluctuation_7A = (energy_7A - mean(energy_7A))/abs(mean(energy_7A));
fluctuation_9A = (energy_9A - mean(energy_9A))/abs(mean(energy_9A));
fluctuation_12A = (energy_12A - mean(energy_12A))/abs(mean(energy_12A));
fluctuation_15A = (energy_15A - mean(energy_15A))/abs(mean(energy_12A));

x = ["7 A", "9 A", "12 A","15 A"];
y = [std(fluctuation_7A),std(fluctuation_9A),std(fluctuation_12A),std(fluctuation_15A)];

figure;
bar(x,y); 
ylabel('总能量相对涨落');
set(gca,'fontsize',10)


