function []=plot_dispersion(Nk,k_norm,name_special_k,nu)
nu=real(nu);
max_nu=max(max(nu));
plot(ones(100,1)*k_norm(1),linspace(0,max_nu*1.1,100),'k-','linewidth',2);
hold on;
for n=1:size(name_special_k,2)-1
    plot(linspace(k_norm(n),k_norm(n+1),Nk),nu(:,(n-1)*Nk+1:n*Nk),'b.');
    plot(ones(100,1)*k_norm(n+1),linspace(0,max_nu*1.1,100),'k-','linewidth',2);
end
set(gca,'xtick',[],'fontsize',12);
ylabel('\nu (THz)','fontsize',12);
axis tight;
for n=1:size(name_special_k,2)
    text(k_norm(n),-max_nu*0.05,name_special_k(n),...
        'interpreter','latex','fontsize',12);
end
