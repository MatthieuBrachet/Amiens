%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 2D CAHN-Hilliard
%
%THE 2d order FINITE DIFFERENCES CASE
%
%The inpainting problem
%
%JPC, 09.03.2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
global initial_type ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUMERICAL DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=63; % number of grid points in each direction 
       %N Must be odd !!!!, ideally 2^p+1
h=1/(N-1);
x=0:h:1;
y=x;
[X,Y]=meshgrid(x,y);
 
%% Time Step size
dt=0.000001;
nper=1;
Nmax=500;

temps=[];

%% Physical data
epsilon=0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ARRAYS
%%%%%%%%%%%%%%%%%%%%%%%%%%
C1=0;
C2=00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construction of the Matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Finrid 
hf=1/(N-1);
DF=gallery('tridiag',N);
DF(1,1)=1;
DF(N,N)=1;
IDF=speye(N,N);
LAPLA=kron(IDF,DF)+kron(DF,IDF);
LAPLA=LAPLA/hf^2;
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time marching scheme matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
m=max(size(LAPLA));

MI11=speye(m,m)+dt*C1*LAPLA ;
MI12=dt*LAPLA;
MI21=-epsilon*LAPLA;
MI22=speye(m,m);
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TIME LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialization
initial_type='triangles';
            
[u0,u0h,lambh]=initial(x,y);
U0=reshape(u0,N*N,1);
U0H=reshape(u0h,N*N,1);
LAMBH=reshape(lambh,N*N,1);
MI11=MI11+dt*spdiags(LAMBH, 0, m, m);
MM=[MI11 MI12; MI21 MI22];
            
U=U0H;
MASS=[mean(U)*h^2];
temps=[0];
k=0;
comp=1;

figure(1)
contourf(X,Y,u0h)
drawnow
pause
while k < Nmax
    clc; disp(abs(k-Nmax)/Nmax*100)

    %Implicit loop in CG
    NLIN=-U.*(1-U.^2)/epsilon;
    %NLIN=U.*(4*U.^2-6*U+2)/epsilon;
    F=[U+dt*C1*LAPLA*U+dt*LAMBH.*(U0);NLIN];
                 
    SOL=MM\F;
    U=SOL(1:m);
    if mod(k,nper)==0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        temps=[temps k*dt];
        MASS=[MASS mean(U)*h^2];

        UFF=reshape(U,N,N);
        contourf(X,Y,UFF)
              
        %film(:,comp)=getframe;
        %fig_placier
        drawnow
        comp=comp+1;
    end  
    k=k+1;
end
figure(2)
plot(temps,MASS)
%post traitement
figure(3)
S=(abs(UFF)>0.5);
                
contourf(X,Y,UFF.*S)
figure(4)
contourf(X,Y,u0h)
figure(5)
contourf(X,Y,S)
%fig_placier
%figure(10)
%movie(film)
%VideoWriter('film.avi');

figure(6)
contourf(X,Y,abs(u0h-S))