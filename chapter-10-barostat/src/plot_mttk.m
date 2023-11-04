clear; close all; font_size=12;
load data.txt;

figure;
histogram(data(:,1),'Normalization','pdf')
xlim([0,10]);
ylim([0,1.5]);
xlabel('x','fontsize',font_size);
ylabel('distribution','fontsize',font_size);
set(gca,'fontsize',font_size);


figure;
histogram(data(:,2),'Normalization','pdf');
xlim([0,10]);
ylim([0,0.5]);
xlabel('L','fontsize',font_size);
ylabel('distribution','fontsize',font_size);
set(gca,'fontsize',font_size);


figure;
histogram(data(:,4),'Normalization','pdf');
xlabel('p','fontsize',font_size);
ylabel('distribution','fontsize',font_size);
set(gca,'fontsize',font_size);


figure;
plot(data(:,3)/data(1,3)-1);
xlabel('time ','fontsize',font_size);
ylabel('E/E_0-1','fontsize',font_size);
set(gca,'fontsize',font_size);
