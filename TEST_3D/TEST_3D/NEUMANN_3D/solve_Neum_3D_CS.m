%
%test
%
clc;clear all;close all;
global K;
%
N=32;
h=1/(N-1);
alpha=1;
beta=1;
%
        x=0:h:1;
        y=x;
        z=x;
        [X,Y,Z]=meshgrid(x,y,z);
%
Nmax=1000;
%
%% 3D Neumann Frequency Matrix
one=ones(N,1);
K=cos(pi*(0:N-1)*h);
%
DX2=kron(one,K);
DY2=DX2';


for j=1:N
     DDX2(:,:,j)=2*(1-DX2);
     DDY2(:,:,j)=2*(1-DY2);
     DDZ2(j,:,:)=2*(1-DX2);
end
      
LAPLA=(DDX2+DDY2+DDZ2)/h^2;



K=alpha+beta*LAPLA;

%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%
%R�solustion du pb de Poisson, sch�mas comapcts
%% solveur classique

[a0,P,Q] = Mlaplacien2(N-2,2);
A=P\Q;
id=speye(N,N);
Ax=kron(kron(A,id),id);
Ay=kron(kron(id,A),id);
Az=kron(kron(id,id),A);
At=sparse(Ax+Ay+Az);
ID=speye(size(At));
MAT=alpha*ID+beta*At;


 %----------------------
        %GMRES Parameters
        %
        max_it=100;
        tol=1.e-12;
        restrt=30;%restarts
        U0=rand(N*N*N,1);
        U33=cos(pi*X).*cos(pi*Y).*cos(pi*Z);
        U3=reshape(U33,N*N*N,1);
        b1=(alpha+3*beta*pi^2)*U3;%ones(n*n*n,1);
        b=reshape(b1,N*N*N,1);
        tic
         [U, error, iter, flag,residual,s] = fouriergmres_neum3D( MAT, U0, b, restrt, max_it, tol );
        tprec=toc;
        
        
        figure(1)
        semilogy(abs(s));
        disp(' *** Error ***** ')
disp([' math. error :' num2str(norm(U-U3,'inf'))])
disp([' time   :' num2str(tprec)])
 disp(' ************** ')