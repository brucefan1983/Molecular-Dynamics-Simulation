clear; close all;
omega=10;dt=0.001:0.001:0.1;
exact_00=cos(omega*dt);
exact_01=sin(omega*dt)/omega;
exact_10=-omega*sin(omega*dt);
exact_11=exact_00;
cayley_prefactor=1./(1+omega*omega*(dt/2).^2);
cayley_00=cayley_prefactor.*(1-omega^2*(dt/2).^2);
cayley_01=cayley_prefactor.*dt;
cayley_10=cayley_prefactor.*(-omega*omega*dt);
cayley_11=cayley_00;

figure;
plot(dt,cayley_00-exact_00,'-','linewidth',2);hold on;
plot(dt,cayley_01-exact_01,'-','linewidth',2);hold on;
legend('00','10');
