clc; clear all; close all;

%% ************************************************************************
global test
tstart=cputime;
%% *** options ************************************************************
film='no';
nper=10;
test=4;
fun='classic'; %% classic or bertozzi

%% *** numerical data
n=40;
[a0,id,AA] = Mlaplacien2(n,2);
Ax=kron(AA,id);
Ay=kron(id,AA);
A=Ax+Ay;
Id=speye(size(A));

[a0,MM,NN] = Mlaplacien2(n,2);
Mx=kron(MM,id);
My=kron(id,MM);
M=Mx+My;

Nx=kron(NN,id);
Ny=kron(id,NN);
N=Nx+Ny;

%% *** space data *********************************************************
h=1/(n+1);
x=[0:h:1]';
[X,Y]=meshgrid(x,x);

%% *** time data **********************************************************
t=0;
ddt=10^-6;
Tmax=0.01;
tau=100;
iter=0;

%% *** mathematical data **************************************************
epsilon=0.05;
[u,unp]=initial_fun(X,Y);
[U]=reshape(u,[],1); [Unp]=reshape(unp,[],1);
U0=U;
W=epsilon*A*U+(1/epsilon).*U.*(U.^2-1);

MAT=[Id tau*ddt*A; -tau*epsilon*A Id];

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
    
    b=[-ddt*LW1-ddt*f; epsilon*LU1-W+(1/epsilon)*NL];
    Z=MAT\b;
    
    U=U+Z(1:(n+2)*(n+2),1);
    W=W+Z((n+2)*(n+2)+1:end,1);
    
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
contourf(X,Y,reshape(triche(U),size(X)))
title('corrected image')
xlabel('x')
ylabel('y')
colorbar

figure(2)
contourf(X,Y,reshape(U0,size(X)))
title('perturbated image')
xlabel('x')
ylabel('y')
colorbar

figure(3)
contourf(X,Y,reshape(U,size(X)))
title('image after Cahn-Hilliard equation')
xlabel('x')
ylabel('y')
colorbar

figure(4)
contourf(X,Y,reshape(triche(U)-Unp,size(X)))
title('error')
xlabel('x')
ylabel('y')
colorbar

