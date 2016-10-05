%
%test
%
clc;clear all;close all;
%
N=32;
h=1/(N-1);
%
        x=0:h:1;
        y=x;
        z=x;
        [X,Y,Z]=meshgrid(x,y,z);
%
epsilon=0.01;
dt=0.01;
tau=1;
Nmax=1000;
Nper=50;
%
expo=exp(-dt/epsilon^2);
expo2=expo^2;
umexpo2=1-expo2;
%
h=1/(N-1);
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

%3D Neumann Laplacian CS matrix
[a0,P,Q] = Mlaplacien2(N-2,4);
A=P\Q;
id=speye(N,N);
Ax=kron(kron(A,id),id);
Ay=kron(kron(id,A),id);
Az=kron(kron(id,id),A);
At=sparse(Ax+Ay+Az);

%
%Iteration matrix for RSS Scheme
%

K=1+tau*dt/2*LAPLA;
EE=speye(N^3,N^3) +dt/2*At;
%Ke=1-dt/2*LAPLA;



%
u=1-2*rand(N*N*N,1);


k=0;
temps=[];
maxu=[];
minu=[];
while k< Nmax

%computation of teh increment by a RRS scheme
auc=reshape(At*u,N,N,N);
auc1=-dt/2*mirt_dctn(auc);
delta=mirt_idctn(auc1./K);
deltastar=reshape(delta,N*N*N,1);
ustar=EE*(u+deltastar);
%second step by exxact integration
u=ustar./sqrt(expo2+umexpo2*ustar.^2);

if  mod(k,Nper)==0
     temps=[temps k*dt];
             figure(1)
             
             timekk=reshape(u,N,N,N);
             show3D
             figure(2)
             isosurface(X,Y,Z,timekk) 
             %contour3(X,Y,Z,U)
             colorbar
             daspect([1,1,1])
             view(3)
             figure(4)
                maxu=[maxu max(u)];
                minu=[minu min(u)];
                plot(temps,maxu,temps,minu)
                xlabel('time')
                legend('max(U)','Min(U)')
             drawnow
             fig_placier ;
end

k=k+1;
end
% u3=reshape(u,N,N,N);
% figure(100)
% iii=(u3>0);
% timek=u3(iii);
% show3D
% figure(200)
% iii=(u3<0);
% timek=u3(iii);
% show3D
