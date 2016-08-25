clc; clear all; close all;
global X Y Z epsilon
% Non-linear RSS

%% time data
ddt=10^-5;
Tmax=1;

%% stabilisation
tau=100;

%% space data
N=40;
h=1/(N+1);
x=[0:h:1]';
[X,Y,Z]=meshgrid(x,x,x);

[a0,MM,A] = Mlaplacien2(N,2);
id=speye(N+2,N+2);
Ax=kron(kron(A,id),id);
Ay=kron(kron(id,A),id);
Az=kron(kron(id,id),A);
Amat=sparse(Ax+Ay+Az);

[a0,MM,NN] = Mlaplacien2(N,4);
id=speye(N+2,N+2);
Mx=kron(kron(MM,id),id);
My=kron(kron(id,MM),id);
Mz=kron(kron(id,id),MM);
Mmat=sparse(Mx+My+Mz);

id=speye(N+2,N+2);
Nx=kron(kron(NN,id),id);
Ny=kron(kron(id,NN),id);
Nz=kron(kron(id,id),NN);
Nmat=sparse(Nx+Ny+Nz);

Id=speye(size(Ax));

%% mathematicals data
t=0;
U=cos(pi*X).*cos(pi*Y).*cos(pi*Z).*exp(sin(3*pi*t));
U=reshape(U,[],1);
epsilon=.5;

%% let's go
ERR=[]; err=0;
TIME=[];
while t<Tmax
    clc; [t err]
    t=t+ddt;
    
     uex=cos(pi*X).*cos(pi*Y).*cos(pi*Z).*exp(sin(3*pi*t));
     Uex=reshape(uex,[],1);
     f=3*pi*cos(3*pi*t)*uex+3*pi^2*uex+(uex.^3-uex)/(epsilon^2);
     f=reshape(f,[],1);
    
     FU=(1/(epsilon^2)).*U.*(U.^2-1);
     DFU=(1/(epsilon^2)).*(3*U.^2-1);
     DF=sparse(diag(DFU));
     
     W=Nmat*U;
     lap=Mmat\W;
     ZZ=(Id+ddt*tau*(Amat+DF))\(-ddt*(lap-FU-f));
     U=ZZ+U;
     
     % erreur
     err=max(abs(Uex-U));
     ERR=[ERR err];
     TIME=[TIME t];
    
end
    
figure(1)
plot(TIME,ERR)
title('NLRSS on Allen-Cahn equation')
xlabel('time')
ylabel('error')