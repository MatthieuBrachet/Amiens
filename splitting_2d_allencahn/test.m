% test : laplacian matrix with neumann BC
% order 4 (hermitian scheme)

clc; clear all; close all;

N=100;
h=1/(N+1);
x=0:h:1;
x=x';

N=N+2;
Qtemp=2*sparse(eye(N,N))-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1);
Qtemp=(6/5)*Qtemp;
Qtemp(1,1)=2681/480;
Qtemp(1,2)=-23/3;
Qtemp(1,3)=113/40;
Qtemp(1,4)=-13/15;
Qtemp(1,5)=59/480;
%N(1,6)=33/40;
Qtemp(end,end)=Qtemp(1,1);
Qtemp(end,end-1)=Qtemp(1,2);
Qtemp(end,end-2)=Qtemp(1,3);
Qtemp(end,end-3)=Qtemp(1,4);
Qtemp(end,end-4)=Qtemp(1,5);
%N(end,end-5)=N(1,6);
Qtemp=(1/h^2)*Qtemp;
Ptemp=sparse(eye(N,N))+1/10*diag(ones(N-1,1),1)+1/10*diag(ones(N-1,1),-1);


u=cos(pi*x);
dduex=pi^2*u;

w=Ptemp\u;
ddu=Qtemp*w;

figure(1)
plot(x,dduex,x,ddu)
legend('exacte','approché')

e=norm(dduex-ddu)





