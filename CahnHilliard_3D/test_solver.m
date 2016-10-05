clc; clear all; close all;
% solver test
%% ************************************************************************
global A Id alpha beta h N

%% *** numerical data
N=32;
[a0,id,AA] = Mlaplacien2(N-2,2);
Ax=kron(kron(AA,id),id);
Ay=kron(kron(id,AA),id);
Az=kron(kron(id,id),AA);
A=Ax+Ay+Az;
Id=speye(size(A));
alpha=1;
beta=10;

MAT=alpha*Id+beta*A^2;

h=1/(N-1);
x=0:h:1;
[X,Y,Z]=meshgrid(x,x,x);


F=rand(size(X));
f=reshape(F,[],1);

SOL2=solv_bilap(F);
sol2=reshape(SOL2,[],1);

sol1=MAT\f;

einf=norm(sol2-sol1,'inf')./norm(sol1,'inf');
