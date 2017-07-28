%% Allen Cahn equation with Euler solver (ADI classic)
clc; clear all; close all; 
epsilon=0.01;
ddt=0.01;
tmax=0.2;
param=20;
tau=2;

n=param-2;
h=1./(n+1);
x=0:h:1;
[X,Y,Z]=meshgrid(x,x,x);
u0=cos(pi*X).*cos(2*pi*Y).*cos(6*pi*Z);


%% *** EULER BLOC
%% space data
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

% initial data
u=reshape(u0,[],1);
time=[]; t=0;

% iterations
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
    
    z=(ID+ddt*A2)\(-ddt*A4u);
    ustar=z+u;
    
    nom=ustar;
    denom=exp(-2*ddt/epsilon^2)+(ustar.^2).*(1-exp(-2*ddt/epsilon^2));
    u=nom./sqrt(denom);

    time=[time t];
end
U1=reshape(u,size(X));
clear u













%% *** RSS BLOC
%% space data
[a0,MM,NN] = Mlaplacien2(n,2);
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

%% initial data
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
U2=reshape(u,size(X));
clear u



%% *** RSS-ADI BLOC
%% space data
[a0,MM,NN] = Mlaplacien2(n,2);
id=speye(size(MM));
Mx=kron(kron(MM,id),id);
Nx=kron(kron(NN,id),id);
My=kron(kron(id,MM),id);
Ny=kron(kron(id,NN),id);
Mz=kron(kron(id,id),MM);
Nz=kron(kron(id,id),NN);
B1=Mx\Nx;
B2=My\Ny;
B3=Mz\Nz;
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
tau=2;

%% initial data
u0=cos(pi*X).*cos(2*pi*Y).*cos(6*Z);
u=reshape(u0,[],1);
time=[]; t=0;

%% iterations
tic;
while t < tmax
    clc; t=t+ddt
    
    wx=Nx*u;
    dux=Mx\wx;
    A4u=dux;
    z=(ID+ddt*tau*B1)\(-ddt*A4u);
    u1=z+u;
    
    wy=Ny*u1;
    duy=My\wy;
    A4u=duy;
    z=(ID+ddt*tau*B2)\(-ddt*A4u);
    u2=z+u1;
    
    wz=Nz*u1;
    duz=Mz\wz;
    A4u=duz;
    z=(ID+ddt*tau*B3)\(-ddt*A4u);
    ustar=z+u2;
    
    nom=ustar;
    denom=exp(-2*ddt/epsilon^2)+(ustar.^2).*(1-exp(-2*ddt/epsilon^2));
    u=nom./sqrt(denom);

    time=[time t];
end
toc
U3=reshape(u,size(X));
clear u



%% *** FIGURES
figure(1)
isosurface(X,Y,Z,U1)
title('Euler')

figure(2)
isosurface(X,Y,Z,U2)
title('RSS')

figure(3)
isosurface(X,Y,Z,U3)
title('RSS-ADI')

fig_placier

disp('error between Euler and RSS:')
max(max(max(abs(U1-U2))))

disp('error between Euler and RSS-ADI:')
max(max(max(abs(U1-U3))))
