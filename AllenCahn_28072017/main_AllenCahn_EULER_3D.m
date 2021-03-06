%% Allen Cahn equation with Euler solver (ADI classic)
clc; clear all; close all; 
epsilon=0.01;
sauvegarde=1;

%% space data
param=16;
n=param-2;
h=1./(n+1);
x=0:h:1;
[X,Y,Z]=meshgrid(x,x,x);

[a0,MM,NN] = Mlaplacien2(n,4);
id=speye(size(MM));
Mx=kron(kron(MM,id),id);
Nx=kron(kron(NN,id),id);
My=kron(kron(id,MM),id);
Ny=kron(kron(id,NN),id);
Mz=kron(kron(id,id),MM);
Nz=kron(kron(id,id),NN);
A2=Mx\Nx+My\Ny+Mz\Nz;
ID=speye(size(A2));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(kron(MM,id),id);
Nx=kron(kron(NN,id),id);
My=kron(kron(id,MM),id);
Ny=kron(kron(id,NN),id);
Mz=kron(kron(id,id),MM);
Nz=kron(kron(id,id),NN);

%% time data and stabilization
tmax=0.2;
ddt=0.01;
tau=1;

%% initial data
u0=cos(pi*X).*cos(2*pi*Y).*cos(6*pi*Z);
u=reshape(u0,[],1);
time=[]; t=0;

%% iterations
tic;
while t < tmax
    clc; t=t+ddt
    
    wx=Nx*u;
    dux=Mx\wx;
    wy=Ny*u;
    duy=My\wy;
    wz=Nz*u;
    duz=My\wz;
    A4u=dux+duy+duz;
    
    z=(ID+ddt*tau*A2)\(-ddt*A4u);
    ustar=z+u;
    
    nom=ustar;
    denom=exp(-2*ddt/epsilon^2)+(ustar.^2).*(1-exp(-2*ddt/epsilon^2));
    u=nom./sqrt(denom);

    time=[time t];
end
toc


U=reshape(u,size(X));
figure(1)
isosurface(X,Y,Z,U)
colorbar

if sauvegarde==1
    print('-dpng', ['surf_N' num2str(param), '_ddt' num2str(1000*ddt), '.png'])
    savefig(['surf_N' num2str(param), '_ddt' num2str(1000*ddt) '.fig']);
end 
