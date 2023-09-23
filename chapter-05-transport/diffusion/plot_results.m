clear all; close all;
load msd.txt;
load vac.txt;

N = size(msd,1);
t = msd(:,1) * 10.18 / 1000; %ps
time_step = t(2) - t(1);
msd = [0; mean(msd(:,2:4),2) / 100]; % nm^2
vac = mean(vac(:,2:4),2) * (1.6/1.66*100); % nm^2/ps^2

D_Einstein = zeros(N, 1);
D_Green_Kubo = zeros(N, 1);

for n = 1 : N
    D_Einstein(n) = (msd(n + 1) - msd(n)) / time_step / 2;
    D_Green_Kubo(n) = sum(vac(1:n)) * time_step;
end
D_Green_Kubo = D_Green_Kubo - 0.5 * time_step * vac(1);

figure(1);
plot(t+time_step, msd(2:end), 'b', 'linewidth', 2);
xlabel('correlation time (ps)');
ylabel('MSD (nm^2)');

figure(2);
plot(t, vac, 'r', 'linewidth', 2);
hold off;
xlabel('correlation time (ps)');
ylabel('VAC (nm^2/ps^2)');

figure(3);
plot(t + time_step, D_Einstein, 'b-', 'linewidth', 4);
hold on;
plot(t + time_step, D_Green_Kubo, 'r--', 'linewidth', 3);
hold off;
xlabel('correlation time (ps)', 'fontsize', 20);
ylabel('D (nm^2/ps)', 'fontsize', 20);
legend('Einstein', 'Green-Kubo');
set(gca, 'fontsize', 20)






