clc; clear all; close all;

%% Allen Cahn 2D extrapol�

tau=3;
eps=1.e-2;

%% espace
N=31;
h=1/(N+1);
x=[0:h:1]';
y=x;
[X,Y]=meshgrid(x,y);


%% temps
t=0;
Tmax=0.1;
dt=1e-5;


%% donnée initiale
%U=reshape(1-2*rand(size(X)),[],1);
U=reshape(cos(pi*X).*cos(pi*Y),[],1);

%% schéma

[a0,A2]=Mlaplacien(N);
[a0,M,N] = Mlaplacien2(N);
I=speye(size(A2));
A2=kron(I,A2)+kron(A2,I);
Mx=kron(I,M);
My=kron(M,I);
Nx=kron(I,N);
Ny=kron(N,I);
Id=speye(size(A2));
global A2 Mx Nx My Ny Id

%% boucle
E=[];T=[];r=1
while t <= Tmax %& r>10^-10
    
    %% rss extrapolé
    
    f=0;
    [V0] = AC1D_iter(U,f,dt/2,tau,eps);
    
    
    f=0;
    [V1] = AC1D_iter(V0,f,dt/2,tau,eps);
    
    
    f=0;
    [V2] = AC1D_iter(U,f,dt,tau,eps);
    
    VV=2*V1-V2;
    
    %% temps
    clc; [t log10(r)]
    t=t+dt;

    %% erreur
    r=norm(VV-U,inf)/dt;
    U=VV;
    
    E=[E r];
    T=[T t];
    
    
    %% figure
    %figure(3)
    %contourf(X,Y,reshape(U,size(X)));
    %colorbar
    
    
    
end
figure(1)
contourf(X,Y,reshape(U,size(X)));
colorbar
xlabel('x')
ylabel('y')
print('-dpng','AC_det1_contour.png')


figure(2)
mesh(X,Y,reshape(U,size(X)));
xlabel('x')
ylabel('y')
print('-dpng','AC_det1_mesh.png')
