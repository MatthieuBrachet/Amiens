clc; clear all; close all;

N=20;
h=1/(N+1);
x=[0:h:1]';
[X,Y,Z]=meshgrid(x,x,x);

U=cos(pi*X).*cos(pi*Y).*cos(pi*Z);
Ubar=reshape(U,[],1);

A=sparse(diag(-2*ones(N+2,1))+diag(ones(N+1,1),1)+diag(ones(N+1,1),-1));
A(1,2)=2; A(end,end-1)=2;
A=-A./(h^2);
Id=speye(N+2,N+2);


Ax=kron(kron(A,Id), Id);
Ay=kron(kron(Id,A), Id);
Az=kron(kron(Id,Id), A);

Lap=Ax+Ay+Az;

ddUbar=Lap*Ubar;
ddUex=3*pi.^2.*Ubar;

E=reshape(abs(ddUbar-ddUex),size(X));

figure(1)
isosurface(X,Y,Z,log10(E))
colorbar


