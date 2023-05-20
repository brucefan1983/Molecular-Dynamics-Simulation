clear;close all;

time_0=[6.515,15.54,33.627,65.971,121.812,1335.68,7194.2];
time_1=[1.425,2.418,3.307,4.957,6.41,22.015,53.86];
time_2=[1.369,2.252,3.612,5.659,8.445,53.63,239.142];
N=[6:10,15,20].^3*4;
figure;
loglog(N,time_0,'d-','linewidth',2); hold on;
plot(N,time_2,'s-','linewidth',2); hold on;
plot(N,time_1,'o-','linewidth',2); hold on;
x=5e3:5e4;
plot(x,1e-6*x.^2,'--','linewidth',2);
plot(x,1e-3*x,':','linewidth',2);
xlabel('Number of atoms');
ylabel('Computation time for 1000 steps');
legend('No neighbor list','O(N^2) algorithm', 'O(N) algorithm)','~N^2','~N');
set(gca,'fontsize',15,'linewidth',1);
axis tight