% ************************************************************************
% ALLEN-CAHN SOLVER
% ************************************************************************
%
% Purpose :    2D numerical simulation to solve the Allen-Cahn equation
%              using spliting between linear and non linear part.
%              the heat equation is solved with the RSS scheme and a space
%              splitting
% Authors :    Brachet Matthieu & Jean-Paul Chehab
% Date    :    2015/07/30
% Last modified : 2015/08/11
%
% ************************************************************************


clc; clear all; close all;

global A A2 Id I P Q tau
global epsilon

%% space data
N=63;
h=1/(N+1);
x=0:h:1;
y=x;
[X,Y]=meshgrid(x,y);

%% time data
dt=1.e-3;
Tmax=1;

%% initial and problem data
tau=1;
epsilon=1.e-2;
U=2*rand(N+2,N+2)-1;
U=reshape(U,[],1);
e=0;

%% Laplacian matrix with Neumann boundaries conditions
A=-2*speye(N+2,N+2)+diag(ones(N+1,1),1)+diag(ones(N+1,1),-1);
A(1,1)=-1; A(end,end)=-1;
Atemp=-sparse(A./h^2);
Id=speye(N+2,N+2);
A2=kron(Atemp,Id)+kron(Id,Atemp);

P=speye(size(A2));
Q=sparse(A2);

I=speye(size(A2));

%% time iterations

for i=1:floor(Tmax/dt)
    
    pause(0.1)
    figure(1)
    contourf(X,Y,reshape(U,N+2,N+2))
    colorbar
    
    clc; [i*dt e]
    
    U1 = nonlinear_part( U, dt/2 );
    U2 = chaleur(U1, dt);
    Unew  = nonlinear_part(U2, dt/2 );
    
    e=norm(Unew-U,inf)./norm(U,inf);
    E(i)=e;
    
    U=Unew;
end


