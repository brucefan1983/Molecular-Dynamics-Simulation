clear; close all;
tic;
E=md(3,4,[4,4,4],5.45*[1,1,1],[1,1,1],1000,1000,20,10,5,80);
toc;

t=(1:50)*0.1; % ps

figure;
plot(t,E(:,2),'^-',t,E(:,1),'v-',t,E(:,3),'o-')
xlabel('Time (ps)','fontsize',15);
ylabel('Energy (eV)','fontsize',15);
legend({'Kinetic','Potential','Total'},'position',[0.4,0.6,0.25,0.1]); 
set(gca,'fontsize',15);

figure;
plot(t,E(:,3)/mean(E(:,3))-1,'o-')
xlabel('Time (ps)','fontsize',15);
ylabel('(E-<E>)/<E>','fontsize',15);
set(gca,'fontsize',15);


