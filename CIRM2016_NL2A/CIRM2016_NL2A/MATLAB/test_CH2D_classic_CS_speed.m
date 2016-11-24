clc; clear all; close all;

  global initial_type ;
%% ************************************************************************
global test
tstart=cputime;
%% *** options ************************************************************
film='no';
nper=10;
test=4;
fun='bertozzi'; %% classic or bertozzi
initial_type='cercles';
%% *** numerical data
n=64;
%% *** space data *********************************************************
h=1/(n+1);
x=[0:h:1];y=x;
[X,Y]=meshgrid(x,x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a0,id,AA] = Mlaplacien2(n,2);
Ax=kron(AA,id);
Ay=kron(id,AA);
A=Ax+Ay;
Id=speye(size(A));

[a0,MM,NN] = Mlaplacien2(n,4);
Mx=kron(MM,id);
My=kron(id,MM);
M=Mx+My;

Nx=kron(NN,id);
Ny=kron(id,NN);
N=Nx+Ny;

A4=M\N; 


%% *** time data **********************************************************
t=0;
ddt=10^-3;
Tmax=0.1;
tau=10;
iter=0;

%% *** mathematical data **************************************************
epsilon=0.05;
%
%INITILISATION
%
%[u,unp]=initial_fun(X,Y);
%[U]=reshape(u,[],1); [Unp]=reshape(unp,[],1);
%U0=U;
%*********************************************
%******************JP CHANGES*****************
m=n+2;
mdim=m*m;
[u0,u0h,lambh]=initial(x,y);
            U0=reshape(u0,m*m,1);
            U0H=reshape(u0h,m*m,1);
            LAMBH=reshape(lambh,m*m,1);
            %MI11=MI11+dt*spdiags(LAMBH, 0, m, m);
            %MM=[MI11 MI12; MI21 MI22];
            %W0H=AA*(U0H+(4*U0H.^2-6*U0H+2)/epsilon);
            U=U0H;
%*********************************************
%*********************************************
W=epsilon*(M\(N*U))+(1/epsilon).*U.*(U.^2-1);

SPL=spdiags(LAMBH,0,m*m,m*m);
%MAT=[(Id+ddt*SPL) tau*ddt*A; -tau*epsilon*A Id];

if strcmp(film,'yes')==1
    nFrames = floor(Tmax/(nper*ddt));
    mov(1:nFrames) = struct('cdata', [],'colormap', []);
    set(gca,'nextplot','replacechildren');
    itermov=0;
end
 
            
%% *** time iterations ****************************************************
while t<=Tmax
    iter=iter+1;
    t=t+ddt;
    clc; disp([iter floor(Tmax/ddt) max(max(U))]);
    
    f=second_membre(X,Y,t,U,U0);
    f=reshape(f,[],1);
    LW1=M\(N*W);
    LU1=M\(N*U);
    
    if strcmp(fun,'classic') == 1
        NL=U.*(U.^2-1);
    elseif strcmp(fun,'bertozzi') == 1
        NL=U.*(4*U.^2-6*U+2);
    else
        error('This non linear part does not exist. Error in ''fun'' charracters. ')
    end
    f=LAMBH.*U0;
    b=[(U+ddt*f);(1/epsilon)*NL];
    
    %
    %Direct solution
    %
    %Z=MAT\b;
    
    %U=U+Z(1:(n+2)*(n+2),1);
    %W=W+Z((n+2)*(n+2)+1:end,1);
    
    %solution via LU bloc
    b1=b(1:mdim);
    b2=b(mdim+1:2*mdim);
    U=(Id+epsilon*ddt*A4^2+ddt*SPL)\(b1-ddt*(A4*b2));
    W=b2+epsilon*A4*U;
    
    
    
    if strcmp(film,'yes')==1 & mod(iter,nper)==0
        figure(1)
        contourf(X,Y,reshape(U,size(X)))
        xlabel('x')
        ylabel('y')
        colorbar
        
        itermov=itermov+1;
        mov(itermov) = getframe(gcf, [0 0 560 420]);
    end
end
disp(['time end : ' num2str(cputime-tstart)])
ref=floor(10000*now);
if strcmp(film,'yes')==1 
    movie2avi(mov, ['ref' num2str(ref) '_date' date '_test' num2str(test) '.avi'], 'compression', 'None');
end

figure(1)
             UUU=reshape(U,m,m); 
                S=(abs(UUU)>0.5);
            contourf(X,Y,S)
            title('Final tresholded solution')
%contourf(X,Y,reshape(triche(U),size(X)))
%title('corrected image')
xlabel('x')
ylabel('y')
colorbar

figure(2)
contourf(X,Y,reshape(U0H,size(X)))
title('Perturbated image - initial situation')
xlabel('x')
ylabel('y')
colorbar

figure(3)
contourf(X,Y,reshape(U,size(X)))
title('Image after Cahn-Hilliard flow')
xlabel('x')
ylabel('y')
colorbar

% figure(4)
% contourf(X,Y,reshape(triche(U)-Unp,size(X)))
% title('error')
% xlabel('x')
% ylabel('y')
% colorbar

