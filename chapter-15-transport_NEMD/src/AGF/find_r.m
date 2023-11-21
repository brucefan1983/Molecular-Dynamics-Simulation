function [r,L,layer_size]=find_r(nxyz)
r0=[1/2,0,0;0,1/6,0;0,1/2,0;1/2,2/3,0];
n0=size(r0,1);
layer_size=nxyz(2)*n0;
N=layer_size*nxyz(1);
a=[1.42*sqrt(3),1.42*3,20];
L=a.*nxyz;
r=zeros(N,3);
n=0;
for nx=0:nxyz(1)-1
    for ny=0:nxyz(2)-1
        for m=1:n0
            n=n+1;
            r(n,:)=a.*([nx,ny,0]+r0(m,:));
        end
    end
end
figure;
plot(r(:,1),r(:,2),'r.','markersize',10);
axis equal;
axis tight;
xlabel('x (A)','fontsize',12);
ylabel('y (A)','fontsize',12);
set(gca,'fontsize',12)
