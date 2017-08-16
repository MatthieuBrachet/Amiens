clc; clear all; close all;
% CODE POUR L'INPAINTING 3D AVEC RSS - 08-15-2017
%% données numériques
N=20;
h=1./(N+1); x=0:h:1; [X,Y,Z]=meshgrid(x,x,x);
ddt=1e-5; t=0;
Tmax=400*ddt; tau=1;
%% données physiques
epsilon=0.01;
lambda=100000; 
test='scroll';
[ u0, u, ind ] = initialdata(X,Y,Z,test);
%% matrices
D=lambda*diag(reshape(ind,[],1));
[a0,id,B] = Mlaplacien2(N,2);
B=kron(kron(B,id),id)+kron(kron(id,B),id)+kron(kron(id,id),B);
ID=speye(size(B));
[a0,M4,N4] = Mlaplacien2(N,4);
Mx=kron(kron(M4,id),id); Nx=kron(kron(N4,id),id);
My=kron(kron(id,M4),id); Ny=kron(kron(id,N4),id);
Mz=kron(kron(id,id),M4); Nz=kron(kron(id,id),N4);
MAT=sparse(ID+D+tau.^2.*ddt*epsilon*B*B);
%% donnees initiales
U0=reshape(u0,[],1); U=reshape(u,[],1);
dx=Nx*U; dxx=Mx\dx;
dy=Ny*U; dyy=My\dy;
dz=Nz*U; dzz=Mz\dz;
Lu=dxx+dyy+dzz;
MU=epsilon*Lu+(1/epsilon).*U.*(1-U.^2);
%% iterations
while t<Tmax
    clc; t=t+ddt
    
    dx=Nx*U; dxx=Mx\dx;
    dy=Ny*U; dyy=My\dy;
    dz=Nz*U; dzz=Mz\dz;
    Lu=dxx+dyy+dzz;
    fu=U.*(1-U.^2);
    
    F1=ddt.*(D*(U0-U)-Lu);
    F2=-MU+epsilon*Lu+(1./epsilon)*fu;
    
    du=MAT\(F1-tau*ddt*B*F2);
    dmu=F2+epsilon*tau*B*du;
    
    U=U+du;
    MU=MU+dmu;
end

u=reshape(U,size(u0));

%% figures
figure(1)
isosurface(X,Y,Z,u0); hold on
colormap gray
isosurface(X,Y,Z,ind); 
colormap bone
hold off
axis([0 1 0 1 0 1])
grid on
title('image initiale avec perturbation')

figure(2)
isosurface(X,Y,Z,u)
colormap gray
axis([0 1 0 1 0 1])
grid on
title('image finale')

ur=sign(u-mean(mean(mean(u))));
figure(3)
isosurface(X,Y,Z,ur)
colormap gray
axis([0 1 0 1 0 1])
grid on
title('image recorrigée')
    
    
fig_placier