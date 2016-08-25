%%%%%%%%%%%%%%%%%%%%%%%%%%
%%test
% FFT Cosinus pour Neumann
%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all


n=62;
h=1/(n+1);
N=n+2;
x=0:h:1;y=x;
[X,Y]=meshgrid(x,y);

ar1=X.*(X-1);
ar2=Y.*(Y-1);
ar=ar1.*ar2;
fa=16;
fapi=fa*pi;
U=cos(fapi*ar);

d1=-fapi*(2*X-1).*ar2.*sin(fapi*ar);
d2=-fapi*(2*Y-1).*ar1.*sin(fapi*ar);
d11=-2*fapi*ar2.*sin(fapi*ar)-fapi^2*(2*X-1).^2.*ar2.^2*cos(fapi*ar);
d22=-2*fapi*ar1.*sin(fapi*ar)-fapi^2*(2*Y-1).^2.*ar1.^2*cos(fapi*ar);

F=U-d11-d22;

k=14;m=5;
U=cos(k*pi*X).*cos(m*pi*Y);
d11=-pi^2*k^2*U;
d22=-pi^2*m^2*U;
F=U-d11-d22;

%
%resolution
%
f=F;%zeros(N,N);
  %spectre de la matrice
  K=1-2*(cos(pi*(0:N-1)'*ones(1,N)*h)+cos(pi*ones(N,1)*(0:N-1)*h)-2)/h^2;
  fhat=dct(dct(f)')';
  uhat=fhat./K;
  u1=idct(idct(uhat)')';
  
  s=max(max(u1))
  max(max(U-u1/s))/h
  mesh(U)
  pause
  mesh(u1)

