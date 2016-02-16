%% test de Navier - Stockes
%% avec fft et extrapolation
clc; clear all; close all;
format short

global cavite

%% *** options ************************************************************
% si sauvegarde = 1 : sauvegarder les graphes;
%               = 0 : ne pas sauvegarder.
sauvegarde = 1;
% si cavite = 1 : cavité raide,
%           = 2 : cavite regularisée.
cavite = 2;
% si film = 1 : faire le film,
%    film = 0 : ne pas faire.
film = 1;

%% données en espace
N=31;
h=1/(N+1);
x=h:h:1-h;
y=x;
[X,Y]=meshgrid(x,y);

%% données en temps
Tmax=1000;
dt=0.026;
t=0;

%% autres données
tau=10;
Re=100;
E=[]

%% fonctions
[Psi,W] = Stokes(N,Re);

%% boucle;
r=1;
while r > 10^-3 & t < Tmax - dt/2

    clc; [t log10(r)]
    
    [bx,by] = bord(Psi);
    [W0] = Convect_diff(dt/2,bx,by,Psi,W,tau,Re);
    [Psi0] = Poisson(N,W0);
    
    [bx,by] = bord(Psi0);
    [W1] = Convect_diff(dt/2,bx,by,Psi0,W0,tau,Re);
    [Psi1] = Poisson(N,W1);
    
    [bx,by] = bord(Psi);
    [W2] = Convect_diff(dt,bx,by,Psi,W,tau,Re);
    [Psi2] = Poisson(N,W2);
    
    PPsi=2*Psi1-Psi2;
    W=2*W1-W2;
    
% PPsi=Psi0;
% W=W0;

    t=t+dt;
    
    r=norm(Psi-PPsi,inf)/dt;
    E=[E r];
    Psi=PPsi;
    
    if log10(r) > 10
        disp('--- ERREUR ---')
        disp('Instabilité')
        disp('Choisir dt plus grand.')
        break 
    end
    
end

%% sauvegarde

save(['NS_N' num2str(N) '_Re' num2str(Re) '_all.mat'])


%% courbes

PP=matrice(Psi);
[IPP1] = level(PP,-1,-0.8*10^-5);
SPP1=PP-IPP1;

hold on
figure(1)
[CC,hh]=contour(X,Y,PP,20);
clabel(CC,hh);
contour(X,Y,SPP1,5);
xlabel('x')
ylabel('y')
title('Psi')
hold off

%print('-dpng',['NS_N' num2str(N) '_Re' num2str(Re) '_1000xdt' num2str(dt*1000) '_psi.png'])

WW=matrice(W);

figure(2)
contour(X,Y,WW,200)
xlabel('x')
ylabel('y')
title('Omega')

%print('-dpng', ['NS_N' num2str(N) '_Re' num2str(Re) '_1000xdt' num2str(dt*1000) '_omega.png'])

%% min et max

[MIN,MAX] = MinMax(PP,X,Y);

disp('max : ')
MAX

disp('min : ')
MIN

%% donn�es

Eextrap=E;
Textrap=[1:length(E)]*dt;

figure(3)
semilogy(Textrap,E)
