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
% si sauvegarde = 1 : sauvegarder les graphes et les données,
%               = 0 : ne pas sauvegarder.
sauvegarde = 1;
% si cavite = 1 : cavité raide,
%           = 2 : cavite regularisée.
cavite = 2;
% si film = 1 : faire le film avec psi,
%    film = 2 : faire le film avec omega,
%    film = 0 : ne pas faire de film.
film = 2;

%% physical data
Re=3200;

%% temps de la comparaison
Tmax=50;
dt=0.0005;
itemax=floor(Tmax./dt);
t=0;

ref=floor(10000*now);

%% *** déclaration des données ********************************************
%% données en espace
N=255;
h=1/(N+1);
x=[h:h:1-h]';
y=x;
[X,Y]=meshgrid(x,y);

%% autres données
tau=1;
E=[];
iter=0;

if film==1 | 2
    % options de film
    nbim=200;%(=5 when dt=0.01)
    
    mkdir(['./video-' date ])
    mov=avifile(['./video-' date '/ref_' num2str(ref) '_NSE.avi'],'compression','None');
    fig=figure;
end

matrice_data
%% *** condition initiale *************************************************
[Psi,W] = Stokes();
%% boucle;
r=1;
while r > 10^-4 & t < Tmax - dt/2
    cc=floor(iter/nbim) - iter/nbim
    clc;  [t cc log10(r)]
    iter=iter+1;
    
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
    
    if film==1
        if floor(iter/nbim) - iter/nbim == 0
           disp('boucle')
           PP=matrice(Psi);
           [IPP1] = level(PP,-1,-0.8*10^-5);
           SPP1=PP-IPP1;

           hold on
           contour(X,Y,PP,20);
           contour(X,Y,SPP1,5);
           xlabel('x')
           ylabel('y')
           title(['Psi - time : ', num2str(t)])
           axis([0 1 0 1])
           
           %film
           frame = getframe(fig, [0 0 560 420]);
           mov = addframe(mov,frame); 
           hold off
           
           close
        end
    elseif film==2
        if floor(iter/nbim) - iter/nbim == 0
           disp('boucle')
           WW=matrice(W);

           hold on
           contour(X,Y,WW,200)
           xlabel('x')
           ylabel('y')
           title(['Omega - time : ', num2str(t)])
           axis([0 1 0 1])
           hold off
           
           %film
           frame = getframe(fig, [0 0 560 420]);
           mov = addframe(mov,frame); 
           hold off
           
           close
        end
    end
    
    
end


%% sauvegarde
if sauvegarde == 1
    mkdir(['./NSE-' date ])
    save(['./NSE-' date '/NS_ref_' num2str(ref), '_all.mat'])
end

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
    print('-dpng',['./NSE-' date '/NS_ref_' num2str(ref), '_psi.png'])
end
WW=matrice(W);

figure(2)
contour(X,Y,WW,200)
xlabel('x')
ylabel('y')
title('Omega')

if sauvegarde == 1 
    print('-dpng',['./NSE-' date '/NS_ref_' num2str(ref), '_omega.png'])
end

%% min et max

[MIN,MAX] = MinMax(PP,X,Y);

disp('max : ')
MAX

disp('min : ')
MIN

if sauvegarde==1
    data = fopen('AAA_NSE_SAVE.txt','a');
    fprintf(data,'%s\n',['date : ', date]);
    fprintf(data,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','---------- numerical data ---------');
    fprintf(data,'%s\n',['grid              : ', num2str(N)] );
    fprintf(data,'%s\n',['time step         : ', num2str(dt)] );
    fprintf(data,'%s\n',['tau               : ', num2str(tau)] );
    fprintf(data,'%s\n','-------- mathematical data --------');
    fprintf(data,'%s\n',['Re                : ', num2str(Re)] );
    fprintf(data,'%s\n',['cavity            : ', num2str(cavite)] );
    fprintf(data,'%s\n','------------ results --------------');
    fprintf(data,'%s\n',['final time        : ', num2str(t)] );
    fprintf(data,'%s\n',['vortex            : '] );
    fprintf(data,'%s\n',['m                 : ', num2str(MIN)] );
    fprintf(data,'%s\n',['m                 : ', num2str(MAX)] );
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','  ');
    fprintf(data,'%s\n','  ');
    fclose(data);
end

if film==1 | 2
    close(fig)
    mov=close(mov);
    
    data = fopen('AAA_VIDEO_SAVE.txt','a');
    fprintf(data,'%s\n',['date : ', date]);
    fprintf(data,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(data,'%s\n',['film : ', num2str(film)]);
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','---------- numerical data ---------');
    fprintf(data,'%s\n',['grid              : ', num2str(N)] );
    fprintf(data,'%s\n',['time step         : ', num2str(dt)] );
    fprintf(data,'%s\n',['tau               : ', num2str(tau)] );
    fprintf(data,'%s\n','-------- mathematical data --------');
    fprintf(data,'%s\n',['Re                : ', num2str(Re)] );
    fprintf(data,'%s\n',['cavity            : ', num2str(cavite)] );
    fprintf(data,'%s\n','------------ results --------------');
    fprintf(data,'%s\n',['final time        : ', num2str(t)] );
    fprintf(data,'%s\n',['vortex            : '] );
    fprintf(data,'%s\n',['m                 : ', num2str(MIN)] );
    fprintf(data,'%s\n',['m                 : ', num2str(MAX)] );
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','  ');
    fprintf(data,'%s\n','  ');
    fclose(data);
end