clc; clear all; close all;

global cavite Re RSS
global sauvegarde
global extrap

%% *** options ************************************************************
% si cavite = 1 : cavité raide,
%           = 2 : cavite regularisée.
cavite = 2;
RSS = 'nonlineaire';
sauvegarde = 1;
extrap = 'no';
%% *** données ************************************************************

N=127;
dt=0.01;
h=1/(N+1);
x=h:h:1-h;
[X,Y]=meshgrid(x,x);
Tmax=10;
Re=3200;

ref=floor(10000*now);
%% *** calculs ************************************************************
tau1=100;
[ T, ~, ~, Ei_pp1, ~, ~, Ei_ww1, PP1, WW1 ] = compNSE(N,dt,Tmax, tau1,Re);


%% graphes

if sauvegarde == 1
    mkdir(['./results-' date ])
    save(['./results-' date '/ref_' num2str(ref) '_data.mat'],'PP1','WW1','T');
end

figure(1)
subplot(121)
semilogy(T,Ei_pp1,'k.')
title(['error on psi - ', num2str(Re)])
legend(['tau = ', num2str(tau1)])
xlabel('time')
ylabel('relative error')

subplot(122)
semilogy(T,Ei_ww1,'k.')
title(['error on omega - ', num2str(Re)])
legend(['tau = ', num2str(tau1)])
xlabel('time')
ylabel('relative error')

if sauvegarde==1
    print('-dpng', ['./results-' date '/ref_' num2str(ref) '_psiomega.png'])
    savefig(['./results-' date '/ref_' num2str(ref) '_psiomega']);

    data = fopen('AAA_RESULTS_SAVE.txt','a');
    fprintf(data,'%s\n',['date : ', date]);
    fprintf(data,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n',['cavity            : ', num2str(cavite)] );
    fprintf(data,'%s\n',['RSS type          : ', RSS]);
    fprintf(data,'%s\n',['extrapolation     : ', extrap]);
    fprintf(data,'%s\n','---------- numerical data ---------');
    fprintf(data,'%s\n',['coarse grid       : ', num2str(N)] );
    fprintf(data,'%s\n',['fine grid         : ', num2str(2*N+1)] );
    fprintf(data,'%s\n',['time step         : ', num2str(dt)] );
    fprintf(data,'%s\n','---------- stabilization ----------');
    fprintf(data,'%s\n',['tau               : ', num2str(tau1)] );
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','  ');
    fprintf(data,'%s\n','  ');
    fclose(data);
end

figure(2)
surf(X,Y,WW1)
title(['omega; tau = ', num2str(tau1)])


if sauvegarde==1
    print('-dpng', ['./results-' date '/ref_' num2str(ref) '_curve_omega.png'])
    savefig(['./results-' date '/ref_' num2str(ref) '_curve_omega']);
end


figure(3)
surf(X,Y,PP1)
title(['psi; tau = ', num2str(tau1)])


if sauvegarde==1
    print('-dpng', ['./results-' date '/ref_' num2str(ref) '_curve_psi.png'])
    savefig(['./results-' date '/ref_' num2str(ref) '_curve_psi']);
end

