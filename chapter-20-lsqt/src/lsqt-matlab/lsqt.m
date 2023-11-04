function [dos,vac,msd,sigma_vac,sigma_msd]=...
    lsqt(Nx,Ny,W,E,E_max,dt,M,flag_vac,flag_msd)
[H,V]=find_H(Nx,Ny,W);
phi=create_state(Nx*Ny);
E=E/E_max;
H=H/E_max;
dos=find_dos(M,E_max,E,H,phi);
vac=zeros(length(dt),length(E));
sigma_vac=zeros(length(dt),length(E));
msd=zeros(length(dt),length(E));
sigma_msd=zeros(length(dt),length(E));
if flag_vac==1
    vac=find_vac(M,E_max,dt*E_max,E,H,V,phi);
    sigma_vac=2*pi*cumtrapz([0;cumsum(dt(1:end-1))],vac);
    for ne=1:length(E)
        vac(:,ne)=vac(:,ne)/dos(ne);
    end
end
if flag_msd==1
    msd=find_msd(M,E_max,dt*E_max,E,H,V/E_max,phi);
    sigma_msd=zeros(length(dt),length(E));
    for ne=1:length(E)
        d_msd=msd(:,ne)-[0;msd(1:end-1,ne)];
        sigma_msd(:,ne)=pi*d_msd./dt;
        msd(:,ne)=msd(:,ne)/dos(ne);
    end
end
