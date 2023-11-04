function energy=get_energy(n_beads,beta)
hbar=1; m=1; lambda=1; dt=0.1; tau_T=100;
omega_n=n_beads/beta/hbar; n_step=1000000; n_step_pimd=1000000;
cayley=true; % cayley is much more stable
C=zeros(n_beads,n_beads);
for j=0:n_beads-1
    for k=0:n_beads-1
        if k==0
            C(j+1,k+1)=sqrt(1/n_beads);
        elseif k<=n_beads/2-1
            C(j+1,k+1)=sqrt(2/n_beads)*cos(2*pi*j*k/n_beads);
        elseif k==n_beads/2
            C(j+1,k+1)=sqrt(1/n_beads)*(-1)^j;
        else
            C(j+1,k+1)=sqrt(2/n_beads)*sin(2*pi*j*k/n_beads);
        end
    end
end
p=linspace(0,0,n_beads); q=linspace(1,1,n_beads); 
energy=zeros(n_step,1);
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
    q_ave=mean(q);
    kinetic_energy=0.5/beta+0.5*m*lambda*lambda*mean((q-q_ave).*q);
    potential_energy=0.5*m*lambda*lambda*mean(q.^2);
    energy(step)=kinetic_energy+potential_energy;
end
energy=mean(energy(end/2+1:end));



