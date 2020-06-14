clear; close all; font_size = 12;
load kappa.txt;

N=200; % number of data for each run
Ns=size(kappa,1)/N % number of independent runs
t=kappa(1:N,1); % time

% running average of the thermal conductivity
kxx=cumsum(reshape(kappa(:,2),N,Ns))./((1:N).'*ones(1,Ns));
kyx=cumsum(reshape(kappa(:,3),N,Ns))./((1:N).'*ones(1,Ns));

disp(['k_xx = (', num2str(mean(kxx(end,:))), ' +- ', ...
    num2str(std(kxx(end,:))/sqrt(Ns)), ') W/mK']);

figure;
plot(t,kxx,'-');
hold on;
plot(t,mean(kxx,2),'--','linewidth',3);
xlabel('Time (ps)','fontsize',font_size);
ylabel('\kappa_{xx} (W/mK)','fontsize',font_size);
title('LJ argon @ T = 20 K from HNEMD');
set(gca,'fontsize',font_size);

figure;
plot(t,kyx,'-');
hold on;
plot(t,mean(kyx,2),'--','linewidth',3);
xlabel('Time (ps)','fontsize',font_size);
ylabel('\kappa_{yx} (W/mK)','fontsize',font_size);
title('LJ argon @ T = 20 K from HNEMD');
set(gca,'fontsize',font_size);

