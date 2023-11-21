function T=find_T1(H00,SigmaL,SigmaR,w,N)
w2=w*w;
H00=w2*eye(N)-H00;
GammaL=imag(SigmaL); 
GammaR=imag(SigmaR);
temp1=(H00-SigmaL-SigmaR)\GammaL;
temp2=(H00-SigmaL'-SigmaR')\GammaR;
T=4*real(trace(temp1*temp2));