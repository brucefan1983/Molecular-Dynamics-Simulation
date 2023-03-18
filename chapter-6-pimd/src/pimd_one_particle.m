clear; close all;
n_beads=8; beta=1; hbar=1; m=1; lambda=5; dt = 0.01; 
omega_n=n_beads/beta/hbar; n_step=200;
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
p=linspace(-2,-1,n_beads); q=linspace(-1,1,n_beads); qq=zeros(n_step,n_beads);
for step=1:n_step
    p=p-(dt/2)*m*lambda*lambda*q; 
    p_normal=p*C; q_normal=q*C;
    for k=0:n_beads-1
        omega_k=2*omega_n*sin(k*pi/n_beads); 
        if k==0
            q_normal(k+1)=(dt/m)*p_normal(k+1)+q_normal(k+1);
        else
            c=cos(omega_k*dt); s=sin(omega_k*dt);
            p_temp=p_normal(k+1);
            q_temp=q_normal(k+1);
            p_normal(k+1)=c*p_temp-m*omega_k*s*q_temp;
            q_normal(k+1)=(1/m/omega_k)*s*p_temp+c*q_temp;
        end
    end
    p=(C*p_normal.').'; q=(C*q_normal.').';
    p=p-(dt/2)*m*lambda*lambda*q;
    qq(step,:)=q;
end
figure;
plot((1:n_step)*dt,mean(qq,2),'r-','linewidth',3);hold on
plot((1:n_step)*dt,qq,'linewidth',0.5);hold on;
xlabel('time');
ylabel('position');
legend('centroid');
set(gca,'fontsize',16);
