clear; close all;
n_beads=128; beta=1; hbar=1; m=1; lambda=1; dt=0.5; tau_T=100;
omega_n=n_beads/beta/hbar; n_step=200000; n_step_pimd=100000;
cayley=true; % cayley is much more stable
C=zeros(n_beads,n_beads);
for j=1:n_beads
    for k=0:n_beads-1
        if k==0
            C(j,k+1)=sqrt(1/n_beads);
        elseif k<=n_beads/2-1
            C(j,k+1)=sqrt(2/n_beads)*cos(2*pi*j*k/n_beads);
        elseif k==n_beads/2
            C(j,k+1)=sqrt(1/n_beads)*(-1)^j;
        else
            C(j,k+1)=sqrt(2/n_beads)*sin(2*pi*j*k/n_beads);
        end
    end
end
p=linspace(0,0,n_beads); q=linspace(1,1,n_beads); 
pp=zeros(n_step,n_beads); qq=zeros(n_step,n_beads);
for step=1:n_step   
    p_normal=p*C;
    c1=exp(-dt*omega_n*sin((0:n_beads-1)*pi/n_beads));
    if step<=n_step_pimd
        c1(1)=exp(-1/2/tau_T);
    end
    c2=sqrt(1-c1.^2);
    p_normal=c1.*p_normal+sqrt(n_beads*m/beta)*c2.*randn(1,n_beads);
    p=(C*p_normal.').';
    p=p-(dt/2)*m*lambda*lambda*q; 
    p_normal=p*C; q_normal=q*C;
    for k=0:n_beads-1
        omega_k=2*omega_n*sin(k*pi/n_beads); 
        if k==0
            q_normal(k+1)=(dt/m)*p_normal(k+1)+q_normal(k+1);
        else
            c=cos(omega_k*dt); s=sin(omega_k*dt);
            if cayley
                c=(1-(omega_k*dt/2)^2)/(1+(omega_k*dt/2)^2);
                s=omega_k*dt/(1+(omega_k*dt/2)^2);
            end
            p_temp=p_normal(k+1);
            q_temp=q_normal(k+1);
            p_normal(k+1)=c*p_temp-m*omega_k*s*q_temp;
            q_normal(k+1)=(1/m/omega_k)*s*p_temp+c*q_temp;
        end
    end
    p=(C*p_normal.').'; q=(C*q_normal.').';
    p=p-(dt/2)*m*lambda*lambda*q;
    p_normal=p*C;
    c1=exp(-dt*omega_n*sin((0:n_beads-1)*pi/n_beads));
    if step<=n_step_pimd
        c1(1)=exp(-1/2/tau_T);
    end
    c2=sqrt(1-c1.^2);
    p_normal=c1.*p_normal+sqrt(n_beads*m/beta)*c2.*randn(1,n_beads);
    p=(C*p_normal.').';
    pp(step,:)=p; qq(step,:)=q;
end
figure;
plot((1:n_step)*dt,qq,'linewidth',0.5);hold on;
plot((1:n_step)*dt,mean(qq,2),'r-','linewidth',3);hold on
xlabel('time');
ylabel('position');
set(gca,'fontsize',16);
figure;
plot(mean(qq(1:end/2,:),2),mean(pp(1:end/2,:),2),'.','markersize',20);
xlabel('position');
ylabel('momentum');
set(gca,'fontsize',16);
figure;
plot(mean(qq(end/2+1:end,:),2),mean(pp(end/2+1:end,:),2),'.','markersize',20);
xlabel('position');
ylabel('momentum');
set(gca,'fontsize',16);


