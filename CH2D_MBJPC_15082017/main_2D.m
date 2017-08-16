clc; clear all; close all;
% CODE POUR L'INPAINTING 2D AVEC RSS - 08-15-2017
%% données numériques
N=63;
h=1./(N+1); x=0:h:1; [X,Y]=meshgrid(x,x);
ddt=0.001; t=0;
Tmax=100*ddt; tau=100;
%% données physiques
epsilon=0.01;
lambda=100000; 
test='triangle';
[ u0, u, ind ] = initialdata(X,Y,test);
%% matrices
D=lambda*diag(reshape(ind,[],1));
[a0,id,B] = Mlaplacien2(N,2);
B=kron(B,id)+kron(id,B);
ID=speye(size(B));
[a0,M4,N4] = Mlaplacien2(N,4);
Mx=kron(id,M4); Nx=kron(id,N4);
My=kron(M4,id); Ny=kron(N4,id);
MAT=ID+D+tau.^2.*ddt*epsilon*B*B;

U0=reshape(u0,[],1); U=reshape(u,[],1);
dx=Nx*U; dxx=Mx\dx;
dy=Ny*U; dyy=My\dy;
Lu=dxx+dyy;
MU=epsilon*Lu+(1/epsilon).*U.*(1-U.^2);
while t<Tmax
    clc; t=t+ddt
    
    dx=Nx*U; dxx=Mx\dx;
    dy=Ny*U; dyy=My\dy;
    Lu=dxx+dyy;
    fu=U.*(1-U.^2);
    
    F1=ddt.*(D*(U0-U)-Lu);
    F2=-MU+epsilon*Lu+(1./epsilon)*fu;
    
    du=MAT\(F1-tau*ddt*B*F2);
    dmu=F2+epsilon*tau*B*du;
    
    U=U+du;
    MU=MU+dmu;
end

u=reshape(U,size(u0));

figure(1)
contourf(X,Y,u)
colormap bone
colorbar
title('image corrigee')

figure(2)
contourf(X,Y,u0)
colormap bone
colorbar
title('image initiale')


ur=(u>0)-(u<0);
figure(3)
contourf(X,Y,ur)
colormap bone
colorbar
title('image recorrigée')
    
    
fig_placier