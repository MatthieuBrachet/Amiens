clc; clear all; close all;
global X Y Z epsilon
%% test Allen-Cahn sans RSS et splitting de Lie

%% time data
ddt=10^-2;
Tmax=1;

%% space data
N=20;
h=1/(N+1);
x=[0:h:1]';
[X,Y,Z]=meshgrid(x,x,x);

[a0,M,A] = Mlaplacien2(N,2);
id=speye(N+2,N+2);

Ax=kron(kron(A,id), id);
Ay=kron(kron(id,A), id);
Az=kron(kron(id,id), A);
Id=speye(size(Ax));

%% mathematicals data
epsilon=0.1;

t=0;
U=sol_exacte(t);
U=reshape(U,[],1);
TIME=[];
ERR=[];
while t<Tmax
    clc; t=t+ddt
    f=sm(t);
    f=reshape(f,[],1);
    
    % step 1
    dt=ddt;
    W=(Id+dt*Ax)\(-dt*Ax*U);
    V1=W+U;
    
    % step 2
    dt=ddt;
    W=(Id+dt*Ay)\(-dt*Ay*V1);
    V2=W+V1;
    
    % step 3
    dt=ddt;
    W=(Id+ddt*Az)\(-ddt*Az*V2+dt*f);
    V3=W+V2;
    
    % step 4
    nom=V3;
    denom=V3.^2+(1-V3.^2).*exp(-2*ddt/(epsilon^2));
    denom=sqrt(denom);
    U=nom./denom;

    %% calcul de l'erreur
    TIME=[TIME t];
    err=max(max(max(abs(reshape(U,size(X))-sol_exacte(t)))));
    ERR=[ERR err];
    
end
figure(1)
isosurface(X,Y,Z,reshape(U,size(X)))

figure(2)
plot(TIME, ERR)

% figure(3)
% subplot(131)
% isosurface(X,Y,Z,reshape(U,size(X)))
% subplot(132)
% isosurface(X,Y,Z,sol_exacte(t))
% subplot(133)
% isosurface(X,Y,Z,reshape(U,size(X))-sol_exacte(t))