clc; clear all; close all;
global X Y Z epsilon
% CNRSS + Splitting de Lie (x-y-z-NL)

%% time data
ddt=10^-4;
Tmax=1;

%% stabilisation
taux=1;
tauy=1;
tauz=1;

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

[a0,MM,NN] = Mlaplacien2(N,2);
id=speye(N+2,N+2);
Mx=kron(kron(MM,id),id);
My=kron(kron(id,MM),id);
Mz=kron(kron(id,id),MM);

id=speye(N+2,N+2);
Nx=kron(kron(NN,id),id);
Ny=kron(kron(id,NN),id);
Nz=kron(kron(id,id),NN);

Zx=Mx\Nx;
Zy=My\Ny;
Zz=Mz\Nz;

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
    
     % step 1
     dt=ddt;
     dU=Nx*U;
     W=(Id+taux*dt/2*Ax)\(-dt/2*(Mx\dU));
%      WW=(Mx+dt/2*Nx)*(U+W);
%      V1=Mx\WW;
    V1=(Id+ddt/2*Zx)*(U+W);
      
     % step 2
     dt=ddt;
     dU=Ny*V1;
     W=(Id+tauy*dt/2*Ay)\(-dt/2*(My\dU));
%      WW=(My+dt/2*Ny)*(V1+W);
%      V2=My\WW;
    V2=(Id+ddt/2*Zy)*(V1+W);
     
     % step 3
     dt=ddt;
     dU=Nz*V2;
     W=(Id+tauz*dt/2*Az)\(-dt/2*(Mz\dU)+dt/2*f);
%      WW=(Mz+dt/2*Nz)*(V2+W);
%      V3=Mz\WW;
    V3=(Id+ddt/2*Zz)*(V2+W);
     
     % step 4
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
title('Linear RSS on Allen-Cahn equation')
xlabel('time')
ylabel('error')