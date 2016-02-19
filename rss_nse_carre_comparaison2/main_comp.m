clc; clear all; close all;

global cavite Re

%% *** options ************************************************************
% si cavite = 1 : cavité raide,
%           = 2 : cavite regularisée.
cavite = 1;

%% *** données ************************************************************

N=127;
dt=0.01;
Tmax=5;
Re=1000;


%% *** calculs ************************************************************
tau1=1;
[ ~, ~, ~, Ei_pp1, ~, ~, Ei_ww1 ] = compNSE(N,dt,Tmax, tau1,Re);

tau2=250;
[ ~, ~, ~, Ei_pp2, ~, ~, Ei_ww2 ] = compNSE(N,dt,Tmax, tau2,Re);

tau3=1000;
[ T, ~, ~, Ei_pp3, ~, ~, Ei_ww3 ] = compNSE(N,dt,Tmax, tau3,Re);

%% graphes
figure(1)
subplot(121)
semilogy(T,Ei_pp1,'k.',T,Ei_pp2,'k--',T,Ei_pp3,'k-')
title(['error on psi - ', num2str(Re)])
legend(['tau = ', num2str(tau1)],['tau = ', num2str(tau2)],['tau = ', num2str(tau3)])
xlabel('time')
ylabel('relative error')

subplot(122)
semilogy(T,Ei_ww1,'k.',T,Ei_ww2,'k--',T,Ei_ww3,'k-')
title(['error on omega - ', num2str(Re)])
legend(['tau = ', num2str(tau1)],['tau = ', num2str(tau2)],['tau = ', num2str(tau3)])
xlabel('time')
ylabel('relative error')

