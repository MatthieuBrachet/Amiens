clc; clear all; close all; 
%% space data
n=50;
h=1./(n+1);
x=0:h:1;
[X,Y]=meshgrid(x,x);

[a0,MM,NN] = Mlaplacien2(n,2);
id=speye(size(MM));
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);
A2=Mx\Nx+My\Ny;
ID=speye(size(A2));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);

%% time data and stabilization
tmax=1;
ddt=0.01;
tau=1;

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
    wy=Ny*u;
    duy=My\wy;
    A4u=dux+duy;
    
    z=(ID+ddt*tau*A2)\(-ddt*A4u+ddt*f);
    u=z+u;
    
    err=[err norm(uex-u,2)./norm(uex,2)];
    time=[time t];
end
toc

figure(1)
semilogy(time,err)