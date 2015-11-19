clc; clear all; close all;

global dt eps
global Mx Nx My Ny A Id
global X Y

%% declaration des variables
% rss data
tau=5;

% time data
Tmax=.1;
dt=0.0001;

% space data
n=200;
h=1/(n+1); x=0:h:1;
[X,Y]=meshgrid(x,x);
[~,I,AA] = Mlaplacien2(n,2);
A=kron(AA,I)+kron(I,AA);
Id=speye(size(A));
[a0,M,N] = Mlaplacien2(n,4);
Mx=kron(M,I);
My=kron(I,M);
Nx=kron(N,I);
Ny=kron(I,N);

% mathematical data
eps=1;
U0=reshape(sol_exacte(0),[],1);


%% resolution d'Allen-Cahn 
[ Uf, T, Eg ] = AC_global(U0, dt, Tmax);
[ Uf, T, Erss ] = AC_RSS(U0, dt, Tmax, tau);

%% figure

figure(1)
semilogy(T,Eg,T,Erss)
legend('classic','residual smoothing scheme',2)