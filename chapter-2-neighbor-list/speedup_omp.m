num=[1,2,4,8,10,12];
t=[26.336,14.213,8.21,7.062,6.371,6.792];
s1=26.336./t;
s2=13.5./t;

figure;
subplot(1,2,1);
plot(num,t,'d-','linewidth',2);
xlabel('Number of threads');
ylabel('Time (s)');
set(gca,'fontsize',15,'linewidth',1);
axis tight;
subplot(1,2,2);
plot(num,s1,'s-','linewidth',2);hold on;
plot(num,s2,'o-','linewidth',2);
xlabel('Number of threads');
ylabel('Speedup');
set(gca,'fontsize',15,'linewidth',1);
legend('相比不用牛顿第三定律的','相比使用牛顿第三定律的')
axis tight;