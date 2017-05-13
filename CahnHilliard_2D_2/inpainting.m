%% Cahn-Hilliard inpainting.
%    start the oct-16-2016
%    authors : Matthieu Brachet & Jean-Paul Chehab
clc; clear all; close all;

global WW XX YY ZZ
global barx bary
global test

test=1;
barx=0.5; bary=0.5;
N=64;
h=1/(N+1);
epsilon=0.1;
lambda=90000;

%% matrices
[a0,id,NN] = Mlaplacien2(N,2);
A2=kron(id,NN)+kron(NN,id);
[a0,MM,NN] = Mlaplacien2(N,4);
A4=MM\NN;
A4=kron(A4,id)+kron(id,A4);
ID=speye(size(A4));

%% time data
ddt=1e-3;
Tmax=1e-2;
t=0;
tau=10000;

%% data
x=0:h:1;
[X,Y]=meshgrid(x,x);
ind=diag(reshape(indicatrice(X,Y,barx,bary),[],1));
[unp,up]=initial_fun(X,Y);
u0=reshape(up,[],1);
u=u0;
w=epsilon*A4*u+(1./epsilon).*u.*(u.^2-1);

%% solver
YY=ID+ddt*lambda*ind;
inD=inv(YY);
ZZ=ddt*tau*A2;
XX=(-epsilon*tau*A2)*inD;
WW=ID-XX*YY;


while t<Tmax
    t=t+ddt;
    clc; disp([t,mean(u)]);
    
    b1=ddt*(-A4*w+lambda*ind*(u0-u));
    b2=epsilon*A4*u-w+(1/epsilon)*u.*(u.^2-1);
    
    [x1,x2]=solver2(b1,b2);
    
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