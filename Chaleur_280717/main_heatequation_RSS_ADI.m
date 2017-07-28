clc; clear all; close all; 
%% space data
param=32;
n=param-2;
h=1./(n+1);
x=0:h:1;
[X,Y]=meshgrid(x,x);

[a0,MM,NN] = Mlaplacien2(n,2);
id=speye(size(MM));
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);
B1=Mx\Nx;
B2=My\Ny;
ID=speye(size(B1));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);

%% time data and stabilization
tmax=1;
ddt=0.001;
tau=2;

%% initial data
u=cos(pi*X).*cos(4*pi*Y);
u=reshape(u,[],1);
err=[]; time=[]; t=0;

%% iterations
tic;
while t < tmax
    clc; t=t+ddt
    uex=cos(pi*X).*cos(4*pi*Y).*exp(sin(pi*t));
    uex=reshape(uex,[],1);
    f=(pi*cos(pi*t)+17*pi^2).*uex;
    
    wx=Nx*u;
    dux=Mx\wx;
    z=(ID+tau*ddt*B1)\(-ddt*dux);
    ustar=z+u;
    
    wy=Ny*ustar;
    duy=My\wy;
    z=(ID+tau*ddt*B2)\(-ddt*duy+ddt*f);
    u=z+ustar;
    
    err=[err norm(uex-u,2)./norm(uex,2)];
    time=[time t];
end
toc

figure(1)
semilogy(time,err)

max(err)