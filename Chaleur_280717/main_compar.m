clc; clear all; close all; 
sauvegarde=1;
%% space data
param=32;
n=param-2;
h=1./(n+1);
x=0:h:1;
[X,Y]=meshgrid(x,x);

%% time data
tmax=1;
ddt=0.0001;





%% *** EULER BLOC
[a0,MM,NN] = Mlaplacien2(n,4);
id=speye(size(MM));
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);
A2=Mx\Nx+My\Ny;
ID=speye(size(A2));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);

%% initial data
u=cos(pi*X).*cos(4*pi*Y);
u=reshape(u,[],1);
err1=[]; time=[]; t=0;

%% iterations
while t < tmax
    clc; t=t+ddt
    uex=cos(pi*X).*cos(4*pi*Y).*exp(sin(pi*t));
    uex=reshape(uex,[],1);
    f=(pi*cos(pi*t)+17*pi^2).*uex;
    
    wx=Nx*u;
    dux=Mx\wx;
    wy=Ny*u;
    duy=My\wy;
    A4u=dux+duy;
    
    z=(ID+ddt*A2)\(-ddt*A4u+ddt*f);
    u=z+u;
    
    err1=[err1 norm(uex-u,2)./norm(uex,2)];
    time=[time t];
end


%% *** RSS BLOC
[a0,MM,NN] = Mlaplacien2(n,2);
id=speye(size(MM));
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);
A2=Mx\Nx+My\Ny;
ID=speye(size(A2));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);

%% initial data
tau=2;
u=cos(pi*X).*cos(4*pi*Y);
u=reshape(u,[],1);
err2=[]; time=[]; t=0;

%% iterations
while t < tmax
    clc; t=t+ddt
    uex=cos(pi*X).*cos(4*pi*Y).*exp(sin(pi*t));
    uex=reshape(uex,[],1);
    f=(pi*cos(pi*t)+17*pi^2).*uex;
    
    wx=Nx*u;
    dux=Mx\wx;
    wy=Ny*u;
    duy=My\wy;
    A4u=dux+duy;
    
    z=(ID+ddt*tau*A2)\(-ddt*A4u+ddt*f);
    u=z+u;
    
    err2=[err2 norm(uex-u,2)./norm(uex,2)];
    time=[time t];
end


%% *** RSS-ADI BLOC
[a0,MM,NN] = Mlaplacien2(n,2);
id=speye(size(MM));
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);
B1=Mx\Nx;
B2=My\Ny;
ID=speye(size(B1));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(MM,id);
Nx=kron(NN,id);
My=kron(id,MM);
Ny=kron(id,NN);

%% initial data
tau=2;
u=cos(pi*X).*cos(4*pi*Y);
u=reshape(u,[],1);
err3=[]; time=[]; t=0;

%% iterations
while t < tmax
    clc; t=t+ddt
    uex=cos(pi*X).*cos(4*pi*Y).*exp(sin(pi*t));
    uex=reshape(uex,[],1);
    f=(pi*cos(pi*t)+17*pi^2).*uex;
    
    wx=Nx*u;
    dux=Mx\wx;
    z=(ID+tau*ddt*B1)\(-ddt*dux);
    ustar=z+u;
    
    wy=Ny*ustar;
    duy=My\wy;
    z=(ID+tau*ddt*B2)\(-ddt*duy+ddt*f);
    u=z+ustar;
    
    err3=[err3 norm(uex-u,2)./norm(uex,2)];
    time=[time t];
end



%% figure
figure(1)
semilogy(time,err1,'k-',time,err2,'k--',time,err3,'k.')
title('2D Heat equation')
xlabel('time')
ylabel('error')
legend('Euler - 4th space order','RSS','RSS-ADI','Location','SouthWest')

if sauvegarde==1
    print('-dpng', ['erreur_N' num2str(param), '_ddt' num2str(1000*ddt), '.png'])
    savefig(['erreur_N' num2str(param), '_ddt' num2str(1000*ddt) '.fig']);
end 