clear; close all;
load temp.txt;
t=(1:size(temp,1))*2/1000;
figure;
plot(t,temp);
xlim([0,40]);
xlabel('Time (ps)','fontsize',12);
ylabel('Temperature (K)','fontsize',12);
set(gca,'fontsize',12);
title('\tau=400 fs');

