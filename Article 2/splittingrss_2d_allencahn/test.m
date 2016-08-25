clc; clear all; close all;

n=1000;
h=1/(n+1);
x=[0:h:1]';
[X,Y]=meshgrid(x,x);
u=reshape(cos(3*pi*X).*cos(3*pi*Y),[],1);

[a0,M,N] = Mlaplacien2(n,4);
Id=speye(size(N));
Mx=kron(M,Id);
My=kron(Id,M);
Nx=kron(N,Id);
Ny=kron(Id,N);

duapp=Mx\(Nx*u)+My\(Ny*u);
ddu=18*pi^2.*u;

DUe=reshape(ddu,size(X));
DUapp=reshape(duapp,size(X));

e=norm(duapp-ddu,inf)

figure(1)
subplot(121)
surf(X,Y,DUapp);
title('approchee')

subplot(122)
surf(X,Y,DUe);
title('exacte')