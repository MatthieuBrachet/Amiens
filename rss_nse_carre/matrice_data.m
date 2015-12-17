%% construction des matrices pour NSE -> main

global N h x Kx Ky bb
global Id Mx My Nx Ny a0
global M1x M1y N1x N1y b0
global cavite


%% *** donnees de bords ***************************************************
b=zeros(N,N);
if cavite == 1
    b(end,:)=1;
elseif cavite == 2
    b(end,:)=(1-(1-2*x).^2).^2;
else
    error('error cavity definition.')
end
bb=vecteur(b);

%% x=0 y=0
A=speye(N,N);
B=zeros(N,N);
B(1,1)=8;
B(1,2)=-3;
B(1,3)=8/9;
B(1,4)=-1/8;

B(end,end)=B(1,1);
B(end,end-1)=B(1,2);
B(end,end-2)=B(1,3);
B(end,end-3)=B(1,4);

Kx=(1/(h^2))*kron(A,B);
Ky=(1/(h^2))*kron(B,A);

%% *** donnes pour les matrices de discrétisation *************************
[a0,Mx,My,Nx,Ny] = Mlaplacien2(N,2);
Id=speye(N,N);

[b0,M1x,M1y,N1x,N1y] = der_prem2(N,2);
