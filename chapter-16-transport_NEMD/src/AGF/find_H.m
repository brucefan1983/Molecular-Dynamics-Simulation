function [H00,H01R,H01L]=find_H(r,r1,r2,L,layer_size,mass,pot)
[NN1,NL1]=find_neighbor(r,1,r1,L);       % for force
[NN2,NL2]=find_neighbor(r,-1,r2,L);      % for harmonic matrix
H00=zeros(layer_size*3,layer_size*3);H01R=zeros(layer_size*3,layer_size*3);
for n1=1:layer_size
    for i2=1:NN2(n1)
        n2=NL2(n1,i2);
        index_1=(n1-1)*3+1:n1*3;
        index_2=(n2-1)*3+1:n2*3;
        H_tmp=find_H1(pot,n1,n2,r,NN1,NL1,L);
        if n2<=layer_size % n1 in layer 0 and n2 in layer 0
            H00(index_1,index_2)=H_tmp;
        elseif n2<=layer_size*2 % n1 in layer 0 and n2 in layer 1
            H01R(index_1,index_2-layer_size*3)=H_tmp;
        end
    end
end
H00=H00/mass;H01R=H01R/mass;H01L=H01R';