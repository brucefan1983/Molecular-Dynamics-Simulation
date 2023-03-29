clear; close all;
m=1; lambda=1; dt=0.01; n_step=1000;
p=0; q=1; 
pp=zeros(n_step,1); qq=zeros(n_step,1);
for step=1:n_step
    p=p-(dt/2)*m*lambda*lambda*q; 
    q=q+(dt/m)*p;
    p=p-(dt/2)*m*lambda*lambda*q;
    pp(step,:)=p; qq(step,:)=q;
end

figure;
plot((1:n_step)*dt,qq,'linewidth',2);hold on;
plot((1:n_step)*dt,pp,'--','linewidth',2);hold on;
xlabel('time');
ylabel('position or momentum');
legend('position','momentum');
set(gca,'fontsize',16);

figure;
plot((1:n_step)*dt,0.5*m*qq.^2,'linewidth',2);hold on;
plot((1:n_step)*dt,pp.^2/m/2,'--','linewidth',2);hold on;
plot((1:n_step)*dt,0.5*m*qq.^2+pp.^2/m/2,'-.','linewidth',3);hold on;
xlabel('time');
ylabel('Energy');
legend('potential','kinetic','total');
set(gca,'fontsize',16);

figure;
plot((1:n_step)*dt,0.5*m*qq.^2+pp.^2/m/2,'-.','linewidth',3);hold on;
xlabel('time');
ylabel('Total Energy');
set(gca,'fontsize',16);

figure;
plot(qq,pp,'.','markersize',20);
xlabel('position');
ylabel('momentum');
set(gca,'fontsize',16);


