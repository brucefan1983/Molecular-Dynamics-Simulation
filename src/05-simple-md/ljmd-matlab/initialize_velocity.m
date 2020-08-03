function [v]=initialize_velocity(K_B,N,D,T,m)
v=rand(N,3)-0.5;
momentum_average=zeros(1,D);
for d=1:D
    momentum_average(d)=sum(v(:,d).*m)/N;
end
for n=1:N
    v(n,:)=v(n,:)-momentum_average/m(n);
end
v=v*sqrt(T*D*K_B*N/sum(m.*sum(v.^2,2))); % scale velocity
