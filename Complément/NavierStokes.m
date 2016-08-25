%% test de Navier - Stockes
%% avec fft
clc; clear all; close all;
global ALAP DX DY;
global a0 ;
global Mx My Nx Ny Id  ;
global RSS ;
%% données en espace
N=127;
h=1/(N+1);
x=h:h:1-h;
y=x;
[X,Y]=meshgrid(x,y);

%% données en temps
Tmax=2000;
dt=0.05;
t=0;

%% autres données
tau=30;
Re=1000;
E=[]

%matrice
matrix_NSE;
[a0,Mx,My,Nx,Ny] = Mlaplacien2(N,2);
Id=speye(size(ALAP));
%type de m�thode (RSS lineaire, RSS non lineaire)
RSS='lineaire';


%% fonctions
[Psi,W] = Stokes(N,Re);

%% boucle;
r=1;
while r > 10^-5 & t < Tmax - dt/2
    clc; [t log10(r)]
    [bx,by] = bord(Psi);
    [W] = Convect_diff(dt,bx,by,Psi,W,tau,Re);
    [Psi2] = Poisson(N,W);
    t=t+dt;
    
   
    r=norm(Psi-Psi2,inf)/dt;
    E=[E r];
    Psi=Psi2;

end


%% courbes

PP=matrice(Psi);
SPP1=level(PP,10^-4,2);
SPP2=level(PP,10^-5,10^-4);

hold on
figure(1)
[CC,hh]=contour(X,Y,PP);
clabel(CC,hh);
[CC,hh]=contour(X,Y,SPP1,2);
clabel(CC,hh);
[CC,hh]=contour(X,Y,SPP2,2);
clabel(CC,hh);
xlabel('x')
ylabel('y')
title('Psi')
hold off

print('-dpng','NS_Re3200_psi_N127.png')

figure(2)
contour(X,Y,matrice(W),200)
xlabel('x')
ylabel('y')
title('Omega')

print('-dpng','NS_Re3200_omega_N127.png')

%% min et max

[MIN,MAX] = MinMax(PP,X,Y);

disp('max : ')
MAX

disp('min : ')
MIN


Re
dt
tau



%% donn�es

Erss=E;
Trss=[1:length(E)]*dt;

figure(3)
semilogy(Trss,Erss)
