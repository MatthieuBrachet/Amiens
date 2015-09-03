% Solve Allen-cahn equation (Neumann BC) on an unit square.
% ---------------------------------------------------------
%
% Strang Splitting is using between non linear and linear part (first step)
% - we use analytic solution for non linear part,
% - we use RSS scheme for heat equation. Preconditionning 
% Laplacian order 4 with Laplacian order 2,
% - we use Strang splitting in heat equation solver.

clc; clear all; close all;
plot_step=0.0001;

global Ax Ay Id I Mx Nx My Ny
global epsilon
global taux tauy

%% space data
N=127;
h=1/(N+1);
x=0:h:1;
y=x;
[X,Y]=meshgrid(x,y);

%% time data
dt=0.001;
Tmax=0.1;
taux=1;
tauy=1;

%% initial and problem data
epsilon=1e-2;
U=reshape(2*rand(N+2,N+2)-1,[],1);
%U=reshape(cos(pi*X).*cos(pi*Y),[],1);

%% Laplacian matrix with Neumann boundaries conditions
A=-2*speye(N+2,N+2)+diag(ones(N+1,1),1)+diag(ones(N+1,1),-1);
A(1,1)=-1; A(end,end)=-1;
Atemp=-sparse(A./h^2);
Id=speye(N+2,N+2);
Ay=kron(Atemp,Id);
Ax=kron(Id,Atemp);

[a0,MM,NN] = Mlaplacien2(N);
Mx=kron(Id,MM);
My=kron(MM,Id);
Nx=kron(Id,NN);
Ny=kron(NN,Id);

I=speye(size(Ax));

%% time iterations

for i=1:floor(Tmax/dt)
    
    %pause(plot_step)
    %figure(1)
    %contourf(X,Y,reshape(U,N+2,N+2))
    %colorbar
    
    clc; i*dt
    
    U1 = nonlinear_part( U, dt/2 );
    U2 = chaleur(U1, dt);
    U  = nonlinear_part(U2, dt/2 );
end

figure(1)
contourf(X,Y,reshape(U,N+2,N+2))
colorbar
    


% the next step is the same with 3D Allen-Cahn solver.