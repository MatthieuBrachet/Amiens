clc; clear all; close all;

%% ************************************************************************
global test

%% *** options ************************************************************
film='no';
test=3;

%% *** numerical data
n=100;
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
Tmax=0.002;
tau=10;

%% *** mathematical data **************************************************
epsilon=0.05;
[u,unp]=initial_fun(X,Y);
[U]=reshape(u,[],1); [Unp]=reshape(unp,[],1);
U0=U;
W=epsilon*A*U+(1/epsilon).*U.*(U.^2-1);

MAT=[Id tau*ddt*A; -tau*epsilon*A Id];

if strcmp(film,'yes')==1
    nFrames = floor(Tmax/ddt);
    mov(1:nFrames) = struct('cdata', [],'colormap', []);
    set(gca,'nextplot','replacechildren');
    iter=1;
end

%% *** time iterations ****************************************************
while t<Tmax
    clc; t=t+ddt
    
    f=second_membre(X,Y,t,U,U0);
    f=reshape(f,[],1);
    LW1=M\(N*W);
    LU1=M\(N*U);
    NL=U.*(U.^2-1);
    b=[-ddt*LW1-ddt*f; epsilon*LU1-W+(1/epsilon)*NL];
    Z=MAT\b;
    
    U=U+Z(1:(n+2)*(n+2),1);
    W=W+Z((n+2)*(n+2)+1:end,1);
    
    if strcmp(film,'yes')==1
        figure(1)
        contourf(X,Y,reshape(U,size(X)))
        xlabel('x')
        ylabel('y')
        colorbar
        
        mov(iter) = getframe(gcf, [0 0 560 420]);
        iter=iter+1;
    end
end
ref=floor(10000*now);
if strcmp(film,'yes')==1
    movie2avi(mov, ['ref' num2str(ref) '_date' date '_test' num2str(test) '.avi'], 'compression', 'None');
end

figure(1)
contourf(X,Y,reshape(triche(U),size(X)))
xlabel('x')
ylabel('y')
colorbar

figure(2)
contourf(X,Y,reshape(U0,size(X)))
xlabel('x')
ylabel('y')
colorbar

figure(3)
contourf(X,Y,reshape(U,size(X)))
xlabel('x')
ylabel('y')
colorbar


