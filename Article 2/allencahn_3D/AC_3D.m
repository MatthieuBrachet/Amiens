clc; clear all; close all;
% Allen Cahn with Lie Splitting and Implicit Euler.
% the solution for the non linear part is explicitly given.
global X Y epsilon

%% time data
ddt=10^-6;
Tmax=0.2;

%% space data
N=20;
h=1/(N+1);
x=[0:h:1]';
[X,Y,Z]=meshgrid(x,x,x);

[a0,MM,A] = Mlaplacien2(N,2);
id=speye(N+2,N+2);

Ax=kron(kron(A,id),id);
Ay=kron(kron(id,A),id);
Az=kron(kron(id,id),A);
Id=speye(size(Ax));

%% mathematicals data
t=0;
U=cos(pi*X).*cos(pi*Y).*cos(pi*Z).*exp(sin(t));
U=reshape(U,[],1);
epsilon=.1;

%% let's go
ERR=[];
TIME=[];
while t<Tmax
    clc; t=t+ddt
    
     uex=cos(pi*X).*cos(pi*Y).*cos(pi*Z).*exp(sin(t));
     Uex=reshape(uex,[],1);
     f=cos(t)*uex+3*pi^2*uex+(uex.^3-uex)/(epsilon^2);
     f=reshape(f,[],1);
    
     % step 1
     dt=ddt;
     W=(Id+dt*Ax)\(-dt*Ax*U);
     V1=W+U;
      
     % step 2
     dt=ddt;
     W=(Id+dt*Ay)\(-dt*Ay*V1);
     V2=W+V1;
     
     % step 3
     dt=ddt;
     W=(Id+dt*Az)\(-dt*Az*V2+dt*f);
     V3=W+V2;
     
     % step 3
     dt=ddt;
     nom=V3;
     denom=V3.^2+(1-V3.^2).*exp(-2*dt/(epsilon^2));
     denom=sqrt(denom);
     U=nom./denom;
     
     % erreur
     err=max(abs(Uex-U));
     ERR=[ERR err];
     TIME=[TIME t];
    
end
    
figure(1)
plot(TIME,ERR)  