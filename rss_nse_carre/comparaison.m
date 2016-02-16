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
% si cavite = 1 : cavité raide,
%           = 2 : cavite regularisée.
cavite = 2;
% si norme = 0 : norme infinie pour le calcul de l'erreur,
%          = p : norme p pour le calcul de l'erreur.
norme = 0;

%% physical data
Re=1000;

%% temps de la comparaison
Tmax=5;
dt=0.0001;
t=0;


%% ************************************************************************
%
%
%                          CALCUL AVEC N=255
%
%
%% ************************************************************************

%% *** déclaration des données ********************************************
%% données en espace
N=63;
h=1/(N+1);
x=[h:h:1-h]';
y=x;
[X,Y]=meshgrid(x,y);

%% autres données
tau=1;
E=[]

matrice_data
%% *** condition initiale *************************************************
[Psi,W] = Stokes();
%% boucle;
r=1;
while r > 10^-5 & t < Tmax - dt/2
    clc; disp('coarse grid :'); [t log10(r)]
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
%% courbes
PP255=matrice(Psi);
WW255=matrice(W);

%% ************************************************************************
%
%
%                          CALCUL AVEC N=511
%
%
%% ************************************************************************

%% *** déclaration des données ********************************************
%% données en espace
N=2*N+1;
h=1/(N+1);
x=[h:h:1-h]';
y=x;
[X,Y]=meshgrid(x,y);

%% autres données
tau=1;
E=[];

matrice_data
%% *** condition initiale *************************************************
[Psi,W] = Stokes();
%% boucle;
r=1;
t=0
while r > 10^-5 & t < Tmax - dt/2
    clc; disp('fine grid :'); [t log10(r)]
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
%% courbes
PP511=matrice(Psi);
WW511=matrice(W);


%% ************************************************************************
%
%
%                          COMPARAISON
%
%
%% ************************************************************************


%% calcul de la restriction
PP511c=PP511(2:2:end-1,2:2:end-1);
WW511c=WW511(2:2:end-1,2:2:end-1);

%% calcul de l'erreur
X=X(2:2:end-1,2:2:end-1);
Y=Y(2:2:end-1,2:2:end-1);

E_PP=abs(PP511c-PP255);
E_WW=abs(WW511c-WW255);

if norme == 0
    e_pp=max(max(E_PP))./max(abs(PP511c));
    e_ww=max(max(E_WW))./max(abs(WW511c));
elseif norme == 2
    e_pp=(sum(sum(abs(PP511c-PP255)).^norme)).^(1/norme)
    e_ww=(sum(sum(abs(WW511c-WW255)).^norme)).^(1/norme)
end

%% sauvegarde
ref=floor(10000*now);
if sauvegarde == 1
    mkdir(['./comp-' date ])
    save(['./comp-' date '/NS_ref_' num2str(ref), '_all.mat'])
end

figure(1)
[CC,hh]=contour(X,Y,E_PP,20);
clabel(CC,hh);
xlabel('x')
ylabel('y')
title('Error on Psi')

if sauvegarde == 1 
    print('-dpng',['./comp-' date '/NS_ref_' num2str(ref), '_erreur_psi.png'])
end

figure(2)
[CC,hh]=contour(X,Y,E_WW,20);
clabel(CC,hh);
xlabel('x')
ylabel('y')
title('Error on Omega')

if sauvegarde == 1 
    print('-dpng',['./comp-' date '/NS_ref_' num2str(ref), '_erreur_omega.png'])
end

figure(3)
surf(X,Y,E_PP);
xlabel('x')
ylabel('y')
title('Error on Psi')

if sauvegarde == 1 
    print('-dpng',['./comp-' date '/NS_ref_' num2str(ref), '_surferreur_psi.png'])
end

figure(4)
surf(X,Y,E_WW);
xlabel('x')
ylabel('y')
title('Error on Omega')

if sauvegarde == 1 
    print('-dpng',['./comp-' date '/NS_ref_' num2str(ref), '_surferreur_omega.png'])
end


if sauvegarde==1
    data = fopen('AAA_COMP_SAVE.txt','a');
    fprintf(data,'%s\n',['date : ', date]);
    fprintf(data,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','---------- numerical data ---------');
    fprintf(data,'%s\n',['fine   grid       : ', num2str(N)] );
    fprintf(data,'%s\n',['coarse grid       : ', num2str((N-1)/2)] );
    fprintf(data,'%s\n',['time step         : ', num2str(dt)] );
    fprintf(data,'%s\n',['tau               : ', num2str(tau)] );
    fprintf(data,'%s\n','-------- mathematical data --------');
    fprintf(data,'%s\n',['Re                : ', num2str(Re)] );
    fprintf(data,'%s\n','------------ results --------------');
    fprintf(data,'%s\n',['final time        : ', num2str(Tmax)] );
    fprintf(data,'%s\n',['norme             : ', num2str(norme)] );
    fprintf(data,'%s\n',['error psi         : ', num2str(e_pp)] );
    fprintf(data,'%s\n',['error omega       : ', num2str(e_ww)] );
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','  ');
    fprintf(data,'%s\n','  ');
    fclose(data);
end