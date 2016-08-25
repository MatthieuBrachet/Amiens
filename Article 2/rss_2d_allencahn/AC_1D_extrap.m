clc; clear all; close all;

%% Allen Cahn 1D extrapolé

tau=5;
eps=10;

%% espace
N=127;
h=1/(N+1);
x=[0:h:1]';
y=x;
[X,Y]=meshgrid(x,y);


%% temps
t=0;
Tmax=5;
dt=0.001;


%% donnée initiale
U=reshape(cos(pi*Y).*cos(pi*X).*exp(sin(0)),[],1);
V=U;


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
E1=[]; E2=[]; T=[];
while t <= Tmax - dt/2
    
    %% rss
    UU=reshape(cos(pi*Y).*cos(pi*X).*exp(sin(t+dt)),[],1);
    f=(cos(t+dt)+2*pi^2).*UU+UU.*(UU.^2+1)/eps;
    [U] = AC1D_iter(U,f,dt,tau,eps);
    
    %% rss extrapolé
    UU=reshape(cos(pi*Y).*cos(pi*X).*exp(sin(t+dt/2)),[],1);
    f=(cos(t+dt/2)+2*pi^2).*UU+UU.*(UU.^2+1)/eps;
    [V0] = AC1D_iter(V,f,dt/2,tau,eps);
    
    UU=reshape(cos(pi*Y).*cos(pi*X).*exp(sin(t+dt)),[],1);
    f=(cos(t+dt)+2*pi^2).*UU+UU.*(UU.^2+1)/eps;
    [V1] = AC1D_iter(V0,f,dt/2,tau,eps);
    
    UU=reshape(cos(pi*Y).*cos(pi*X).*exp(sin(t+dt)),[],1);
    f=(cos(t+dt)+2*pi^2).*UU+UU.*(UU.^2+1)/eps;
    [V2] = AC1D_iter(V,f,dt,tau,eps);
    
    V=2*V1-V2;
    
    %% temps
    clc; t=t+dt

    %% erreur
    E1=[E1 norm(UU-U,inf)];
    E2=[E2 norm(UU-V,inf)];
    T=[T t];
    
    
end

figure(1)
semilogy(T,E1,T,E2)
legend('rss','extrap')
xlabel('temps')
ylabel('erreur')

print('-dpng','AC_N127_tau5_eps10.png')