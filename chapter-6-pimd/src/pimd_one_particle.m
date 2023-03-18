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
    p_normal=zeros(1,n_beads);
    q_normal=zeros(1,n_beads);
    for k=1:n_beads
        for j=1:n_beads
            p_normal(k)=p_normal(k)+p(j)*C(j,k);
            q_normal(k)=q_normal(k)+q(j)*C(j,k);
        end
    end
    
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
    for j=1:n_beads
        p(j)=0;q(j)=0;
        for k=1:n_beads
            p(j)=p(j)+p_normal(k)*C(j,k);
            q(j)=q(j)+q_normal(k)*C(j,k);
        end
    end
    
    p=p-(dt/2)*m*lambda*lambda*q;
    qq(step,:)=q;
end
figure;
plot((1:n_step)*dt,qq,'linewidth',1);hold on;
plot((1:n_step)*dt,mean(qq,2),'--','linewidth',3);
xlabel('time');
ylabel('position');
legend('bead-0','bead-1','bead-2','bead-3','bead-4','bead-5','bead-6','bead-7','centroid');
set(gca,'fontsize',16);
