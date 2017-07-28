%% Allen Cahn equation with Euler solver (ADI classic)
clc; clear all; close all; 
epsilon=0.01;
ddt=0.01;
tmax=0.2;
param=32;
tau=2;
sauvegarde=1;


%% *** EULER BLOC
%% space data
n=param-2;
h=1./(n+1);
x=0:h:1;
[X,Y]=meshgrid(x,x);

[a0,MM,NN] = Mlaplacien2(n,4);
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

%% initial data
u0=cos(pi*X).*cos(2*pi*Y);
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
    A4u=dux+duy;
    
    z=(ID+ddt*A2)\(-ddt*A4u);
    ustar=z+u;
    
    nom=ustar;
    denom=exp(-2*ddt/epsilon^2)+(ustar.^2).*(1-exp(-2*ddt/epsilon^2));
    u=nom./sqrt(denom);

    time=[time t];
end
toc
U1=reshape(u,size(X));


%% *** RSS BLOC
%% space data
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
A2=Mx\Nx+My\Ny;
ID=speye(size(A2));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);

%% initial data
u=cos(pi*X).*cos(2*pi*Y);
u=reshape(u,[],1);
time=[]; t=0;

%% iterations
while t < tmax
    clc; t=t+ddt
    
    wx=Nx*u;
    dux=Mx\wx;
    wy=Ny*u;
    duy=My\wy;
    A4u=dux+duy;
    
    z=(ID+ddt*tau*A2)\(-ddt*A4u);
    ustar=z+u;
    
    nom=ustar;
    denom=exp(-2*ddt/epsilon^2)+(ustar.^2).*(1-exp(-2*ddt/epsilon^2));
    u=nom./sqrt(denom);
    
    time=[time t];
end
U2=reshape(u,size(X));



%% *** RSS-ADI BLOC
%% space data
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
ID=speye(size(B2));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);


%% initial data
u=cos(pi*X).*cos(2*pi*Y);
u=reshape(u,[],1);
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
    ustar=z+u1;
    
    nom=ustar;
    denom=exp(-2*ddt/epsilon^2)+(ustar.^2).*(1-exp(-2*ddt/epsilon^2));
    u=nom./sqrt(denom);

    time=[time t];
end
toc
U3=reshape(u,size(X));



%% *** FIGURES
figure(1)
contourf(X,Y,U1)
title('Euler')
colorbar

figure(2)
contourf(X,Y,U2)
title('RSS')
colorbar

figure(3)
contourf(X,Y,U3)
title('RSS-ADI')
colorbar

figure(4)
contourf(X,Y,abs(U1-U2))
title('RSS against Euler')
colorbar

figure(5)
contourf(X,Y,abs(U1-U3))
title('RSS-ADI against Euler')
colorbar

figure(6)
contourf(X,Y,abs(U2-U3))
title('RSS against RSS-ADI')
colorbar

fig_placier