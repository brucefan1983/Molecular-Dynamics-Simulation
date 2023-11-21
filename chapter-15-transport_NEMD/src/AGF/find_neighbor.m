function [NN,NL]=find_neighbor(r,r1,r2,L)
N=size(r,1);
NN=zeros(N,1);
NL=zeros(N,100);
for n1=1:N
    for n2=1:N
        r12=r(n2,:)-r(n1,:);
        r12=r12-round(r12./L).*L;
        d12=norm(r12);
        if (d12>r1)&&(d12<r2)
            NN(n1)=NN(n1)+1;
            NL(n1,NN(n1))=n2;
        end
    end
end
NL=NL(:,1:max(NN));