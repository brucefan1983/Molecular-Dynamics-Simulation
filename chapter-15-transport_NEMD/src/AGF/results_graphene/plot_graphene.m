clear; close all
load T; load nu;
temp=300;
A=3.35*1.42*3*24*1.0e-20;%m^2
kB=1.38e-23;
g=1.0e3*kB*T/A;
h=6.63e-34;
x=h*nu*1.0e12/kB/temp;
fx=x.^2.*exp(x)./(exp(x)-1).^2;
g_quantum=g.*fx;


figure;
plot(nu,g,'linewidth',2);
hold on;
plot(nu,g_quantum,'--','linewidth',2);
xlabel('\omega/2\pi (THz)','fontsize',12);
ylabel('G(\omega) (GW/m^2/K/THz)','fontsize',12);
xlim([0,55]);
ylim([0,0.5]);
set(gca,'fontsize',12,'ticklength',get(gca,'ticklength')*2);
legend('Classical','Quantum');
text(2,0.46,'Pristine graphene @ 300 K','fontsize',12)

trapz(nu,g)
trapz(nu,g_quantum)
