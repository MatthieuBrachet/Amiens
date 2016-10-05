clc; clear all; close all;
% solver test
%% ************************************************************************
global A Id alpha beta h N LAPLA

N=32;
h=1/(N-1);
x=0:h:1;
[X,Y,Z]=meshgrid(x,x,x);


[a0,id,AA] = Mlaplacien2(N-2,2);
Ax=kron(kron(AA,id),id);
Ay=kron(kron(id,AA),id);
Az=kron(kron(id,id),AA);
A=Ax+Ay+Az;
Id=speye(size(A));
alpha=1;
beta=10;

MAT=alpha*Id-beta*A;

DX2=cos(pi*(0:N-1)'*ones(1,N)*h);
DY2=cos(pi*ones(N,1)*(0:N-1)*h);
for j=1:N
    DDX2(:,:,j)=2*(1-DX2);
    DDY2(:,:,j)=2*(1-DY2);
end
for k=1:N
    DDZ2(k,:,:)=2*(1-DX2);
end
      
LAPLA=(DDX2+DDY2+DDZ2)/h^2;

      
F=rand(size(X));
f=reshape(F,[],1);

SOL2=solv_poiss_neumann_3D(F);
sol2=reshape(SOL2,[],1);

sol1=MAT\f;

einf=norm(sol2-sol1,'inf')./norm(sol1,'inf');