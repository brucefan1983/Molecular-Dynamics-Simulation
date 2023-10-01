function [r,basis_label,type,L]=find_r(r0,basis_label0,type0,nxyz,axyz)
L=axyz.*nxyz;
M=size(r0,1);
N=prod(nxyz)*M;
count=0;
r=zeros(N,3);
basis_label=zeros(N,1);
type=zeros(N,1);
for i1=1:nxyz(1)
    for i2=1:nxyz(2)
        for i3=1:nxyz(3)
            for index=1:M
                x=axyz(1)*(i1-1)+r0(index,1);
                y=axyz(2)*(i2-1)+r0(index,2);
                z=axyz(3)*(i3-1)+r0(index,3);
                count=count+1;
                r(count,:)=[x,y,z];
                basis_label(count)=basis_label0(index);
                type(count)=type0(index);
            end
        end
    end
end
