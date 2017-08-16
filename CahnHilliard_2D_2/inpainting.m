%% Cahn-Hilliard inpainting : 2D
%    start the oct-16-2016
%    authors : Matthieu Brachet & Jean-Paul Chehab
clc; clear all; close all;

global XX YY ZZ
global barx bary
global test

test=1;
barx=0.5; bary=0.1;
N=64;
h=1/(N+1);
epsilon=5;
lambda=1000000;

%% matrices
[a0,id,NN] = Mlaplacien2(N,2);
A2=kron(id,NN)+kron(NN,id);
[a0,MM,NN] = Mlaplacien2(N,2);
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);
ID=speye(size(Mx));

%% time data
ddt=1e-3;
Tmax=0.2;
t=0;
tau=10000;

%% data
x=0:h:1;
[X,Y]=meshgrid(x,x);
ind=diag(reshape(indicatrice(X,Y,barx,bary),[],1));
[unp,up]=initial_fun(X,Y);
u0=reshape(up,[],1);
u=u0;
ux=Nx*u; ux=Mx\ux; 
uy=Ny*u; uy=My\uy; 
w=epsilon*(ux+uy)+(1./epsilon).*u.*(u.^2-1);

%% solver
XX=sparse(ID+lambda*ind+tau*tau*ddt*epsilon*A2^2);
YY=-tau*ddt*A2;
ZZ=epsilon*tau*A2;


while t<Tmax
    t=t+ddt;
    clc; disp([t,mean(u)]);
    
    wx=Nx*w; wx=Mx\wx; 
    wy=Ny*w; wy=My\wy; 
    b1=ddt*(-wx-wy+lambda*ind*(u0-u));
    
    ux=Nx*u; ux=Mx\ux; 
    uy=Ny*u; uy=My\uy; 
    b2=epsilon*(ux+uy)-w+(1/epsilon)*u.*(u.^2-1);
    
    [x1,x2]=solver3(b1,b2);
    
    u=u+x1;
    w=w+x2;
end

im_init=reshape(unp,size(X));
im_pert=up;
im_finale=reshape(u,size(X));

figure(1);
contourf(X,Y,im_init);
title('initial image')
colorbar

figure(2);
contourf(X,Y,im_pert);
title('pertubated image')
colorbar

figure(3)
contourf(X,Y,im_finale);
title('corrected image')
colorbar

figure(4)
contourf(X,Y,triche(im_finale));
title('corrected and adjust image')
colorbar

fig_placier