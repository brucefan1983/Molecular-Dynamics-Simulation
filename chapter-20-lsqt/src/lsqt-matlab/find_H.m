function [H,V]=find_H(Nx,Ny,W)
N=Nx*Ny;
row=zeros(N*4,1);
col=zeros(N*4,1);
Hij=-ones(N*4,1);
Vij=zeros(N*4,1);
count=0;
for nx=0:Nx-1
    for ny=0:Ny-1
        index_center=nx*Ny+ny+1;
        if nx>0
            index_left=(nx-1)*Ny+ny+1;
            count=count+1;
            row(count)=index_center;
            col(count)=index_left;
            Vij(count)=1i;
        end
        if nx<Nx-1
            index_right=(nx+1)*Ny+ny+1;
            count=count+1;
            row(count)=index_center;
            col(count)=index_right;
            Vij(count)=-1i;
        end
        if ny<Ny-1
            index_up=nx*Ny+ny+2;
            count=count+1;
            row(count)=index_center;
            col(count)=index_up;
        end
        if ny>0
            index_down=nx*Ny+ny;
            count=count+1;
            row(count)=index_center;
            col(count)=index_down;
        end
    end
end
row=row(1:count);
col=col(1:count);
Hij=Hij(1:count);
Vij=Vij(1:count);
H=sparse(row,col,Hij,N,N,count);
U=sparse(1:N,1:N,(rand(N,1)-0.5)*W,N,N,N);
H=H+U;
V=sparse(row,col,Vij,N,N,count);
end
