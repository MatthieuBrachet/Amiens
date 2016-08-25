clc; clear all; close all;
global X Y epsilon

%% time data
ddt=10^-4;
Tmax=0.5;

%% space data
N=50;
h=1/(N+1);
x=[0:h:1]';
[X,Y]=meshgrid(x,x);

A=sparse(diag(-2*ones(N+2,1))+diag(ones(N+1,1),1)+diag(ones(N+1,1),-1));
A(1,2)=2; A(end,end-1)=2;
A=-A./(h^2);
id=speye(N+2,N+2);

Ax=kron(A,id);
Ay=kron(id,A);
Id=speye(size(Ax));

%% mathematicals data
t=0;
U=cos(pi*X).*cos(pi*Y).*exp(sin(t));
U=reshape(U,[],1);
epsilon=.1;

%% let's go
ERR=[];
TIME=[];
while t<Tmax
    clc; t=t+ddt
    
     uex=cos(pi*X).*cos(pi*Y).*exp(sin(t));
     Uex=reshape(uex,[],1);
     f=cos(t)*uex+2*pi^2*uex+(uex.^3-uex)/(epsilon^2);
     f=reshape(f,[],1);
    
     % step 1
     dt=ddt;
     Z=(Id+dt*Ax)\(-dt*Ax*U);
     V1=Z+U;
      
     % step 2
     dt=ddt;
     Z=(Id+dt*Ay)\(-dt*Ay*V1+dt*f);
     V2=Z+V1;
     
     % step 3
     dt=ddt;
     nom=V2;
     denom=V2.^2+(1-V2.^2).*exp(-2*dt/(epsilon^2));
     denom=sqrt(denom);
     U=nom./denom;
     
     % erreur
     err=max(abs(Uex-U));
     ERR=[ERR err];
     TIME=[TIME t];
    
end
    
figure(1)
plot(TIME,ERR)  