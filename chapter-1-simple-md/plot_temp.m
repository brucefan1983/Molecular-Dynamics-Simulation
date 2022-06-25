clear; close all
load energy.txt;
timeStep = 5/1000; % ps
sampleInterval = 100;
timeInterval = timeStep * sampleInterval;
numData = size(energy, 1);
time = (1 : numData) * timeInterval;
figure;

plot(time, energy(:,1), '-', 'linewidth', 2);hold on;
xlabel('Time (ps)', 'fontsize', 15);
ylabel('Temperature (K)', 'fontsize', 15);
set(gca, 'fontsize', 15);

