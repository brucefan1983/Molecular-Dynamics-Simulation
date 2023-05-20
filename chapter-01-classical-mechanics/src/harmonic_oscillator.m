clear; close all;
m=1; k=1; dt=0.01; n_step=1000;
v=0; x=1; 
v_vector=zeros(n_step,1); x_vector=zeros(n_step,1);
for step=1:n_step
    v=v-(dt/2)*k*x; 
    x=x+dt*v;
    v=v-(dt/2)*k*x; 
    v_vector(step,:)=v; x_vector(step,:)=x;
end

figure;
plot((1:n_step)*dt,x_vector,'linewidth',2);hold on;
plot((1:n_step)*dt,v_vector,'--','linewidth',2);hold on;
xlabel('time');
ylabel('position or velocity');
legend('position','velocity');
set(gca,'fontsize',16);

potential=0.5*k*x_vector.^2;
kinetic=0.5*m*v_vector.^2;
figure;
plot((1:n_step)*dt,potential,'linewidth',2);hold on;
plot((1:n_step)*dt,kinetic,'--','linewidth',2);hold on;
plot((1:n_step)*dt,potential+kinetic,'-.','linewidth',3);hold on;
xlabel('time');
ylabel('Energy');
legend('potential','kinetic','total');
set(gca,'fontsize',16);

figure;
plot((1:n_step)*dt,potential+kinetic,'-.','linewidth',3);hold on;
xlabel('time');
ylabel('Total Energy');
set(gca,'fontsize',16);

figure;
plot(x_vector,v_vector,'.','markersize',20);
xlabel('position');
ylabel('momentum');
set(gca,'fontsize',16);


