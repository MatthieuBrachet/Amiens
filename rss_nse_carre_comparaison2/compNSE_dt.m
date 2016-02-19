clc; clear all; close all;

global cavite Re

%% *** options ************************************************************
% si cavite = 1 : cavité raide,
%           = 2 : cavite regularisée.
cavite = 1;

%% *** données ************************************************************

N=31;
tau=1; 
Tmax=5;
Re=100;


%% *** calculs ************************************************************
dt1=0.016;
[ T1, ~, ~, Ei_pp1, ~, ~, Ei_ww1 ] = compNSE(N,dt1,Tmax,tau,Re);

dt2=0.01;
[ T2, ~, ~, Ei_pp2, ~, ~, Ei_ww2 ] = compNSE(N,dt2,Tmax,tau,Re);

dt3=0.001;
[ T3, ~, ~, Ei_pp3, ~, ~, Ei_ww3 ] = compNSE(N,dt3,Tmax,tau,Re);

%% graphes
figure(1)
subplot(121)
title(['error on psi - ', num2str(Re)])
semilogy(T1,Ei_pp1,'k.',T2,Ei_pp2,'k--',T3,Ei_pp3,'k-')
legend(['dt = ', num2str(dt1)],['dt = ', num2str(dt2)],['dt = ', num2str(dt3)])
xlabel('time')
ylabel('relative error')

subplot(122)
title(['error on omega - ', num2str(Re)])
semilogy(T1,Ei_ww1,'k.',T2,Ei_ww2,'k--',T3,Ei_ww3,'k-')
legend(['dt = ', num2str(dt1)],['dt = ', num2str(dt2)],['dt = ', num2str(dt3)])
xlabel('time')
ylabel('relative error')