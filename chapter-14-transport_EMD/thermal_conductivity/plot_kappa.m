clear;close all; font_size=15;
load kappa.txt; % Data from the C code

% number of correlation data for one simulation; 
% should be consistent with the calculations using the C code
M = 200; 

% I usually use 10-100 simulations to get high-quality data
number_of_simulations = size(kappa, 1) / M 

t = kappa(1 : M, 1); % correlation time

% For isotropic 3D systems, k = (kx + ky + kz) / 3
hac = reshape(mean(kappa(:, 2:4), 2), M, number_of_simulations);
rtc = reshape(mean(kappa(:, 5:7), 2), M, number_of_simulations);

% ensemble average
hac_ave = mean(hac, 2);
rtc_ave = mean(rtc, 2);

% Normalized heat current auto-correlation function
figure
for n = 1 : number_of_simulations
    plot(t, hac(:, n) / hac(1, n), 'color', 0.5 * [1 1 1]);
    hold on;
end
plot(t, hac_ave ./ hac_ave(1), 'linewidth', 3);
xlabel('Correlation Time (ps)', 'fontsize', font_size);
ylabel('HCACF (Normalized)', 'fontsize', font_size);
set(gca,'fontsize', font_size);

% Running thermal conductivity
figure
for n = 1 : number_of_simulations
    plot(t, rtc, 'color', 0.5 * [1 1 1]);
    hold on;
end
plot(t, rtc_ave, 'linewidth', 3);
xlabel('Correlation Time (ps)', 'fontsize', font_size);
ylabel('\kappa (W/mK)', 'fontsize', font_size);
set(gca,'fontsize', font_size);

% Average over an appropriate time block after visual inspection
kappa_converged = mean(rtc(M/2+1 : M, :));

% Report the mean value and an error estimate (standard error)
kappa_average = mean(kappa_converged)
kappa_error = std(kappa_converged) / sqrt(number_of_simulations)