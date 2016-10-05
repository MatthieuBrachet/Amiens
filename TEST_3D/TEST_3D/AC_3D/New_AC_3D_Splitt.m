%
%test
%
clc;clear all;close all;
%
N=16;
h=1/(N-1);
%
        x=0:h:1;
        y=x;
        z=x;
        [X,Y,Z]=meshgrid(x,y,z);
%
epsilon=0.01;
dt=0.1;
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



K=1+dt/2*LAPLA;
Ke=1-dt/2*LAPLA;
KK=1+dt*LAPLA;


%
u=1-2*rand(N,N,N);
uc=mirt_dctn(u);

k=0;
temps=[];
maxu=[];
minu=[];
while k< Nmax


ustar=mirt_idctn((Ke./K).*uc);

u=ustar./sqrt(expo2+umexpo2*ustar.^2);

if  mod(k,Nper)==0
     temps=[temps k*dt];
             figure(1)
             timekk=u;
             show3D
             figure(2)
             isosurface(X,Y,Z,u) 
             %contour3(X,Y,Z,U)
             colorbar
             daspect([1,1,1])
             view(3)
             figure(4)
                maxu=[maxu max(max(max(u)))];
                minu=[minu min(min(min(u)))];
                plot(temps,maxu,temps,minu)
                xlabel('time')
                legend('max(U)','Min(U)')
             drawnow
             fig_placier ;
end
unext=u;
uc=mirt_dctn(u);
k=k+1;
end
figure(100)
iii=(u>0);
v=u.*iii;
timekk=1-v;
show3D
figure(200)
iii=(u<0);
v=u.*iii;
timekk=1-v;
show3D
