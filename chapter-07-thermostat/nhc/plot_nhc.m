clear; close all; font_size=20;
load data.txt;

figure;
plot(data(:,1),data(:,2),'.');
axis equal;
xlabel('x (arbitrary unit)','fontsize',font_size);
ylabel('p (arbitrary unit)','fontsize',font_size);
xlim([-5,5]);
ylim([-5,5]);

figure;
hist(data(:,1),31);
xlabel('x (arbitrary unit)','fontsize',font_size);
ylabel('count','fontsize',font_size);
xlim([-5,5]);

figure;
hist(data(:,2),31);
xlabel('p (arbitrary unit)','fontsize',font_size);
ylabel('count','fontsize',font_size);
xlim([-5,5]);
