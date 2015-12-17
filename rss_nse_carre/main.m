clear all; close all; clc;

%% ************************************************************************
% Résolution de l'équationd e Navier-Stokes (Vorticité, Fct de Courant)
% par le schéma RSS en différences finies aux ordre 2 et 4 sur un carré.
%
% *** authors
%        Matthieu Brachet
%        Jean-Paul Chehab
%
%% ************************************************************************

%% *** variables globales *************************************************
global N h X Y x y
global Re
global cavite

%% *** options ************************************************************
% si sauvegarde = 1 : sauvegarder les graphes;
%               = 0 : ne pas sauvegarder.
sauvegarde = 1;
% si cavite = 1 : cavité raide
%           = 2 : cavite regularisée
cavite = 1;
%% *** déclaration des données ********************************************
%% données en espace
N=63;
h=1/(N+1);
x=[h:h:1-h]';
y=x;
[X,Y]=meshgrid(x,y);

%% données en temps
Tmax=1000;
dt=0.001;
t=0;

%% autres données
tau=1;
Re=1000;
E=[]

matrice_data
%% *** condition initiale *************************************************
[Psi,W] = Stokes();
%% boucle;
r=1;
while r > 10^-5 & t < Tmax - dt/2
    clc; [t log10(r)]
    [bx,by] = bord(Psi);
    [W0] = Convect_diff(dt/2,bx,by,Psi,W,tau);
    [Psi0] = Poisson(N,W0);
    
    [bx,by] = bord(Psi0);
    [W1] = Convect_diff(dt/2,bx,by,Psi0,W0,tau);
    [Psi1] = Poisson(N,W1);
    
    [bx,by] = bord(Psi);
    [W2] = Convect_diff(dt,bx,by,Psi,W,tau);
    [Psi2] = Poisson(N,W2);
    
    PPsi=2*Psi1-Psi2;
    W=2*W1-W2;
    
    t=t+dt;
    r=norm(Psi-PPsi,inf)/dt;
    E=[E r];
    Psi=PPsi;
    if log10(r) > 10
        error('Instability : dt, the time step must be smaller or \tau bigger.')
        break 
    end
end
%% sauvegarde
ref=floor(10000*now);
if sauvegarde == 1
    mkdir(['./results-' date ])
    save(['./results-' date '/NS_ref_' num2str(ref), '_all.mat'])
end

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

if sauvegarde == 1 
    print('-dpng',['./results-' date '/NS_ref_' num2str(ref), '_psi.png'])
end

WW=matrice(W);
figure(2)
contour(X,Y,WW,200)
xlabel('x')
ylabel('y')
title('Omega')
if sauvegarde == 1
    print('-dpng', ['./results-' date '/NS_ref_' num2str(ref), '_omega.png'])
end

%% min et max
[MIN,MAX] = MinMax(PP,X,Y);
disp('max : ')
MAX
disp('min : ')
MIN
%% analysed e convergence vers un etat stationnaire.

Eextrap=E;
Textrap=[1:length(E)]*dt;

figure(3)
semilogy(Textrap,E)
if sauvegarde == 1
    print('-dpng', ['./results-' date '/NS_ref_' num2str(ref), '_convergence.png'])
end

if sauvegarde==1
    data = fopen('AAA_RESULTS_SAVE.txt','a');
    fprintf(data,'%s\n',['date : ', date]);
    fprintf(data,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','---------- numerical data ---------');
    fprintf(data,'%s\n',['number of points  : ', num2str(N)] );
    fprintf(data,'%s\n',['time step         : ', num2str(dt)] );
    fprintf(data,'%s\n',['tau               : ', num2str(tau)] );
    fprintf(data,'%s\n','-------- mathematical data --------');
    fprintf(data,'%s\n',['Re                : ', num2str(Re)] );
    fprintf(data,'%s\n','------------ results --------------');
    fprintf(data,'%s\n',['final time        : ', num2str(Textrap(end))] );
    fprintf(data,'%s\n',['convergence data  : ', num2str(E(end))] );
    fprintf(data,'%s\n',['primary vortex    : ', num2str(MAX)] );
    fprintf(data,'%s\n',['secondary vortex  : ', num2str(MIN)] );
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','  ');
    fprintf(data,'%s\n','  ');
    fclose(data);
end