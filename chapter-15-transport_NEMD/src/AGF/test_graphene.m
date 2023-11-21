clear;close all;font_size=12;line_width=2; 
nu=0.2:0.2:52; % THz
w=nu*2*pi*10.18/1000;
r1=2;r2=2.6;mass=12;
[r0,L,layer_size]=find_r([3,24,1]);
pot.model=0; % only tersoff now
pot.para=[1393.6,430,3.5333,2.2407,1.5724e-7,0.72751,38049,4.3484,-0.93,1.8,2.1];
T=find_T(w,r0,r1,r2,L,layer_size,mass,pot);
figure;   
plot(nu,T,'r-','linewidth',line_width);
xlabel('$\omega$ (THz)','fontsize',font_size,'interpreter','latex');
ylabel('Transmission','fontsize',font_size);
xlim([0,nu(end)*1.1]);
ylim([-0.1,max(T)*1.1]);
set(gca,'fontsize',font_size)

save('T','T');
save('nu','nu');