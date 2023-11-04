function [SigmaL,SigmaR]=find_Sigma(H00,H01R,H01L,w)
gsR=find_g00(H00,H01R,w);
gsL=find_g00(H00,H01L,w);
SigmaR=H01R*gsR*H01R';
SigmaL=H01L*gsL*H01L';