clc; clear all; close all;

%% test chaleur

%% time data
t=0;
ddt=0.0001;
tau=1;
Tmax=1;

%% space data and functions
NN=100;
h=1/(NN+1);
x=[0:h:1]';

u=cos(pi*x)*exp(-t);
uex=u;

%% matricial data
[a0,M,N] = Mlaplacien2(NN,4);
A4=M\N;
[a0,M,N] = Mlaplacien2(NN,2);
A2=M\N;
Id=speye(size(M));
MAT=Id - (ddt/2) * inv(Id+0.5*tau*ddt*A2)*A4;
%% time step
E=[0];
while t<Tmax & E(end)<10^10
    t=t+ddt;
    clc; disp([t E(end)]);
    
    
    f1=-uex+pi^2*uex;
    uex=cos(pi*x)*exp(-t);
    f2=-uex+pi^2*uex;
    f=0.5*(f1+f2);
    
    % Euler Implicite
    %unew=(Id+ddt*A2)\(ddt*f+u);
    
    % CN
    %unew=(Id+ddt/2*A2)\(ddt*f+(Id-ddt/2*A2)*u);
    
    % CNRSS1
    w=ddt*f+(Id-0.5*ddt*A4)*u;
    unew=MAT*w;
    
    % CNRSS2
    %delta=(Id+0.5*tau*ddt*A2)\(-ddt/2*A4*u);
    %unew=(Id+ddt/2*A4)*(u+delta);
    
    e=norm(uex-unew);
    E=[E e];
    
    u=unew;
end

figure(1)
plot(E)