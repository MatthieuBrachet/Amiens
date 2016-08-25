clc; clear all; close all;
global X Y

%% time data
ddt=10^-4;
Tmax=0.4;

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

t=0;
U=cos(pi*X).*cos(pi*Y).*exp(sin(t));
U=reshape(U,[],1);

ERR=[];
TIME=[];
while t<Tmax
    clc; t=t+ddt
    
    Uex=cos(pi*X).*cos(pi*Y).*exp(sin(t));
    Uex=reshape(Uex,[],1);
    f=cos(t)*Uex+2*pi^2*Uex;
    
     % step 1
     dt=ddt;
     Z=(Id+dt*Ax)\(-dt*Ax*U);
     V1=Z+U;
      
     % step 2
     dt=ddt;
     Z=(Id+dt*Ay)\(-dt*Ay*V1+dt*f);
     U=Z+V1;
    
     % erreur
     err=max(abs(Uex-U));
     ERR=[ERR err];
     TIME=[TIME t];
    
end
    
figure(1)
plot(TIME,ERR)  