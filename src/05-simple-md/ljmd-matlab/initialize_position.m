function [r]=initialize_position(N,D,n0,nxyz,a) % FCC crystal
r0=[0,0,0;0,0.5,0.5;0.5,0,0.5;0.5,0.5,0];r=zeros(N,D);n=0;
for nx=0:nxyz(1)-1
    for ny=0:nxyz(2)-1
        for nz=0:nxyz(3)-1
            for m=1:n0
                n=n+1;r(n,:)=a.*([nx,ny,nz]+r0(m,:));
            end
        end
    end
end
