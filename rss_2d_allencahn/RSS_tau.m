clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mesure de l'erreur dans AC par RSS extrapol� ou non
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% donn�s

eps=0.5;

%% espace
N=63;
h=1/(N+1);
x=[0:h:1]';
y=x;
[X,Y]=meshgrid(x,y);

%% temps
Tmax=1;
dt=0.0001;

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


clc; 1

%% Allen Cahn 2D extrapol�

%% donnée initiale
t=0;
U=reshape(cos(pi*X).*cos(pi*Y).*exp(sin(t)),[],1);

%% boucle
E=[];T=[];r=1
tau=1;
while t <= Tmax - dt/2 & r>10^-10
    
    %% rss extrapolé
    
    F=f(X,Y,t+dt/2,eps);
    [V0] = AC1D_iter(U,F,dt/2,tau,eps);
    
    
    F=f(X,Y,t+dt,eps);
    [V1] = AC1D_iter(V0,F,dt/2,tau,eps);
    
    
    F=f(X,Y,t+dt,eps);
    [V2] = AC1D_iter(U,F,dt,tau,eps);
    
    U=2*V1-V2;
    
    %% temps
    %clc; [t log10(r)]
    t=t+dt;

    %% erreur
    Uex=reshape(cos(pi*X).*cos(pi*Y).*exp(sin(t)),[],1);
    r=norm(U-Uex,inf);
    
    E=[E r];
    T=[T t];
end



2



%% Allen Cahn 2D extrapol�

%% donnée initiale
t=0;
U=reshape(cos(pi*X).*cos(pi*Y).*exp(sin(t)),[],1);

%% boucle
E1=[];T1=[];r=1
tau=3;
while t <= Tmax - dt/2 & r>10^-10
    
    %% rss extrapolé
    
    F=f(X,Y,t+dt/2,eps);
    [V0] = AC1D_iter(U,F,dt/2,tau,eps);
    
    
    F=f(X,Y,t+dt,eps);
    [V1] = AC1D_iter(V0,F,dt/2,tau,eps);
    
    
    F=f(X,Y,t+dt,eps);
    [V2] = AC1D_iter(U,F,dt,tau,eps);
    
    U=2*V1-V2;
    
    %% temps
    %clc; [t log10(r)]
    t=t+dt;

    %% erreur
    Uex=reshape(cos(pi*X).*cos(pi*Y).*exp(sin(t)),[],1);
    r=norm(U-Uex,inf);
    
    E1=[E1 r];
    T1=[T1 t];
end

3



%% Allen Cahn 2D

%% donnée initiale
t=0;
U=reshape(cos(pi*X).*cos(pi*Y).*exp(sin(t)),[],1);

%% boucle
E2=[];T2=[];r=1
tau=10;
while t <= Tmax - dt/2 & r>10^-10
    
    %% rss extrapolé
    
    F=f(X,Y,t+dt/2,eps);
    [V0] = AC1D_iter(U,F,dt/2,tau,eps);
    
    
    F=f(X,Y,t+dt,eps);
    [V1] = AC1D_iter(V0,F,dt/2,tau,eps);
    
    
    F=f(X,Y,t+dt,eps);
    [V2] = AC1D_iter(U,F,dt,tau,eps);
    
    U=2*V1-V2;
    
    %% temps
    %clc; [t log10(r)]
    t=t+dt;

    %% erreur
    Uex=reshape(cos(pi*X).*cos(pi*Y).*exp(sin(t)),[],1);
    r=norm(U-Uex,inf);
    
    E2=[E2 r];
    T2=[T2 t];
end




4


%% Allen Cahn 2D

%% donnée initiale
t=0;
U=reshape(cos(pi*X).*cos(pi*Y).*exp(sin(t)),[],1);

%% boucle
E3=[];T3=[];r=1
tau=100;
while t <= Tmax - dt/2 & r>10^-10
    
    %% rss extrapolé
    
    F=f(X,Y,t+dt/2,eps);
    [V0] = AC1D_iter(U,F,dt/2,tau,eps);
    
    
    F=f(X,Y,t+dt,eps);
    [V1] = AC1D_iter(V0,F,dt/2,tau,eps);
    
    
    F=f(X,Y,t+dt,eps);
    [V2] = AC1D_iter(U,F,dt,tau,eps);
    
    U=2*V1-V2;
    
    %% temps
    %clc; [t log10(r)]
    t=t+dt;

    %% erreur
    Uex=reshape(cos(pi*X).*cos(pi*Y).*exp(sin(t)),[],1);
    r=norm(U-Uex,inf);
    
    E3=[E3 r];
    T3=[T3 t];
end







figure(1)
semilogy(T,E,'b--',T1,E1,'b',T2,E2,'r',T3,E3,'g');
legend('RSS tau=1','RSS tau=3','RSS tau=10','RSS tau=100')
xlabel('time')
ylabel('error')
print('-dpng','AC_tau1.png')