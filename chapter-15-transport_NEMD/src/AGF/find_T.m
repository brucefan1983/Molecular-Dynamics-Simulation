function [T]=find_T(w,r,r1,r2,L,layer_size,mass,pot)
Nw=length(w);
T=zeros(1,Nw);  
[H00,H01R,H01L]=find_H(r,r1,r2,L,layer_size,mass,pot);
for n=1:Nw
    w_n=w(n);
    [SigmaL,SigmaR]=find_Sigma(H00,H01R,H01L,w_n);
    T(n)=find_T1(H00,SigmaL,SigmaR,w_n,layer_size*3);
end