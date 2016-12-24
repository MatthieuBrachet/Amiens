%% Cahn-Hilliard inpainting 3D.
%    start the nov-28-2016
%    authors : Matthieu Brachet & Jean-Paul Chehab
clc; clear all; close all;

global WW XX YY ZZ
global barx bary barz centre
global test

%% video options
film='yes';

%%
test=3;
barx=.8; bary=.1; barz=.8;
centre=[.5 .5 .5];
N=20;
h=1/(N+1);
epsilon=0.05;
lambda=100000;

%% matrix
[a0,id,NN] = Mlaplacien2(N,2);
A2=sparse(kron(kron(NN,id),id)+kron(kron(id,NN),id)+kron(kron(id,id),NN));
[a0,MM,NN] = Mlaplacien2(N,4);
A4=MM\NN;
A4=sparse(kron(kron(A4,id),id)+kron(kron(id,A4),id)+kron(kron(id,id),A4));
ID=speye(size(A4));

%% time data
ddt=1e-7;
Tmax=100*ddt;
t=0;
tau=1;

%% data
x=0:h:1;
[X,Y,Z]=meshgrid(x,x,x);
ind=diag(reshape(indicatrice(X,Y,Z,barx,bary,barz),[],1));
[unp,up]=initial_fun(X,Y,Z);
u0=reshape(up,[],1);
u=u0;
w=epsilon*A4*u+(1./epsilon).*u.*(u.^2-1);

%% solver
YY=sparse(ID+ddt*lambda*ind);
inD=inv(YY);
ZZ=ddt*tau*A2;
XX=(-epsilon*tau*A2)*inD;
WW=ID-XX*YY;

iter=0;
while t<Tmax
    t=t+ddt;
    iter=iter+1;
    iter=floor(t/ddt);
    clc; disp([iter,mean(u)]);
    
    b1=ddt*(-A4*w+lambda*ind*(u0-u));
    b2=epsilon*A4*u-w+(1/epsilon)*u.*(u.^2-1);
    
    [x1,x2]=solver2(b1,b2);
    
    u=u+x1;
    w=w+x2;
    
    if strcmp(film,'yes')==1
        im_v(:,:,:,iter)=reshape(u,size(X));
    end
    
end


im_init=reshape(unp,size(X));
im_pert=up;
im_finale=reshape(u,size(X));
im_finale01=reshape(triche(u),size(X));

figure(1)
isosurface(X,Y,Z,im_init);
title('initial function')
axis([0 1 0 1 0 1])
xlabel('x')
ylabel('y')
zlabel('z')
grid on

figure(2)
isosurface(X,Y,Z,im_pert);
title('perturbated image')
axis([0 1 0 1 0 1])
xlabel('x')
ylabel('y')
zlabel('z')
grid on

figure(3)
isosurface(X,Y,Z,im_finale)
title('corrected image')
axis([0 1 0 1 0 1])
xlabel('x')
ylabel('y')
zlabel('z')
grid on

figure(4)
isosurface(X,Y,Z,im_finale01)
title('re-corrected image')
axis([0 1 0 1 0 1])
xlabel('x')
ylabel('y')
zlabel('z')
grid on

fig_placier