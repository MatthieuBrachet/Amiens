%% Allen Cahn equation with RSS solver (ADI classic)
%error('ne fonctionne pas encore')
clc; clear all; close all; 
epsilon=0.01;

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
A2=Mx\Nx+My\Ny;
ID=speye(size(A2));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);

%% time data and stabilization
tmax=0.2;
ddt=0.0001;
tau=1;

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
    wy=Ny*u;
    duy=My\wy;
    A4u=dux+duy;
    
    DF=0*(-3*u.^2+1)/(epsilon^2);
    DF=sparse(diag(DF));
    F=(u-u.^3)./epsilon.^2;
   
    z=(ID+ddt*tau*A2+ddt*DF)\(-ddt*A4u-ddt*F);
    u=z+u;

    time=[time t];
end
toc


U=reshape(u,size(X));
figure(1)
contourf(X,Y,U)
colorbar


