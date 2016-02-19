function [ T, E1_pp, E2_pp, Ei_pp, E1_ww, E2_ww, Ei_ww ] = compNSE(N,dt,Tmax, tau,Re)

%% données en espace 1
N1=N;
%h1=1/(N1+1);
%x1=h1:h1:1-h1;
%y1=x1;
%[X1,Y1]=meshgrid(x1,y1);

%% données en espace 2
N2=2*N1+1;
%h2=1/(N2+1);
%x2=h2:h2:1-h2;
%y2=x2;
%[X2,Y2]=meshgrid(x2,y2);

%% données en temps
t=0;


%% historique
Ei_pp=[];
Ei_ww=[];
E2_pp=[];
E2_ww=[];
E1_pp=[];
E1_ww=[];
T=[];

%% fonctions
[Psi1,W1] = Stokes(N1,Re);
[Psi2,W2] = Stokes(N2,Re);

%% boucle;
while t < Tmax
    clc; [tau t]
    
    %% discretisation d'espace N1
    [bx1,by1] = bord(Psi1);
    [W01] = Convect_diff(dt/2,bx1,by1,Psi1,W1,tau,Re);
    [Psi01] = Poisson(N1,W01);
    
    [bx1,by1] = bord(Psi01);
    [W11] = Convect_diff(dt/2,bx1,by1,Psi01,W01,tau,Re);
    [Psi11] = Poisson(N1,W11);
    
    [bx1,by1] = bord(Psi1);
    [W21] = Convect_diff(dt,bx1,by1,Psi1,W1,tau,Re);
    [Psi21] = Poisson(N1,W21);
    
    Psi1=2*Psi11-Psi21;
    W1=2*W11-W21;
    
    %% discretisation d'espace N2
    [bx2,by2] = bord(Psi2);
    [W02] = Convect_diff(dt/2,bx2,by2,Psi2,W2,tau,Re);
    [Psi02] = Poisson(N2,W02);
    
    [bx2,by2] = bord(Psi02);
    [W12] = Convect_diff(dt/2,bx2,by2,Psi02,W02,tau,Re);
    [Psi12] = Poisson(N2,W12);
    
    [bx2,by2] = bord(Psi2);
    [W22] = Convect_diff(dt,bx2,by2,Psi2,W2,tau,Re);
    [Psi22] = Poisson(N2,W22);
    
    Psi2=2*Psi12-Psi22;
    W2=2*W12-W22;
    
    %% calcul de l'erreur
    PP1=matrice(Psi1);
    WW1=matrice(W1);
    PP2=matrice(Psi2);
    WW2=matrice(W2);
    
    PP2c=PP2(2:2:end-1,2:2:end-1);
    WW2c=WW2(2:2:end-1,2:2:end-1);

    E_PP=abs(PP2c-PP1);
    E_WW=abs(WW2c-WW1);

    % norme infinie
    e_pp=max(max(E_PP))./max(max(abs(PP2c)));
    e_ww=max(max(E_WW))./max(max(abs(WW2c)));
        
    Ei_pp=[Ei_pp e_pp];
    Ei_ww=[Ei_ww e_ww];
    
    if Ei_pp(end) > 100 || Ei_ww(end) > 100
        error('instability : choose different dt or tau')
        break
    end
        
    % norme 2
    norme = 2;
    e_pp=(sum(sum(abs(PP2c-PP1).^norme))).^(1/norme);
    e_pp=e_pp./((sum(sum(abs(PP2c).^norme))).^(1/norme));
    e_ww=(sum(sum(abs(WW2c-WW1).^norme))).^(1/norme);
    e_ww=e_ww./((sum(sum(abs(WW2c).^norme))).^(1/norme));
        
    E2_pp=[E2_pp e_pp];
    E2_ww=[E2_ww e_ww];
        
    % norme 1
    norme = 1;
    e_pp=(sum(sum(abs(PP2c-PP1).^norme))).^(1/norme);
    e_pp=e_pp./((sum(sum(abs(PP2c).^norme))).^(1/norme));
    e_ww=(sum(sum(abs(WW2c-WW1).^norme))).^(1/norme);
    e_ww=e_ww./((sum(sum(abs(WW2c).^norme))).^(1/norme));
        
    E1_pp=[E1_pp e_pp];
    E1_ww=[E1_ww e_ww];
    
    %% temps
    T=[T t];
    t=t+dt;
end


end

