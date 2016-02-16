function [FF] = NL(Psi,W,bx,by)
%% partie non linéaire de navier_stokes
N=sqrt(length(Psi));
[b0,M1x,M1y,N1x,N1y] = der_prem2(N,2);

aa0=zeros(N,N);
aa0(:,1)=b0;
aa0(:,end)=-b0;
a0y=vecteur(aa0');

aa0=zeros(N,N);
aa0(:,1)=b0;
aa0(:,end)=-b0;
a0x=vecteur(aa0);

dwx=M1x\(N1x*W+a0x.*bx);
dwy=M1y\(N1y*W+a0y.*by);

dpsix=M1x\(N1x*Psi);
dpsiy=M1y\(N1y*Psi);

FF=dpsiy.*dwx-dpsix.*dwy;


end