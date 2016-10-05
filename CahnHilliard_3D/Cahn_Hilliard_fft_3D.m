clc; clear all; close all;
% with 3D solver (fft).
%% ************************************************************************
global test
global A Id nb h
tstart=cputime;
%% *** options ************************************************************
film='no';
nper=5;
test=2;
fun='classic'; %% classic or bertozzi

%% *** numerical data
nb=32;
[~,id,AA] = Mlaplacien2(nb-2,2);
Ax=kron(kron(AA,id),id);
Ay=kron(kron(id,AA),id);
Az=kron(kron(id,id),AA);
A=Ax+Ay+Az;
Id=speye(size(A));

[a0,MM,NN] = Mlaplacien2(nb-2,4);
Mx=kron(kron(MM,id),id);
My=kron(kron(id,MM),id);
Mz=kron(kron(id,id),MM);
M=sparse(Mx+My+Mz);

Nx=kron(kron(NN,id),id);
Ny=kron(kron(id,NN),id);
Nz=kron(kron(id,id),NN);
N=sparse(Nx+Ny+Nz);

%% *** space data *********************************************************
h=1/(nb-1);
x=0:h:1;
[X,Y,Z]=meshgrid(x,x,x);

%% *** time data **********************************************************
t=0;
ddt=10^-7;
Tmax=1000*ddt;
tau=1;
iter=0;

%% *** mathematical data **************************************************
epsilon=0.05;

[u,unp]=initial_fun(X,Y,Z);
[U]=reshape(u,[],1); [Unp]=reshape(unp,[],1);
U0=U;
W=epsilon*A*U+(1/epsilon).*U.*(U.^2-1);


if strcmp(film,'yes')==1
    nFrames = floor(Tmax/(nper*ddt));
    mov(1:nFrames) = struct('cdata', [],'colormap', []);
    set(gca,'nextplot','replacechildren');
    itermov=0;
    
    figure(1)
    isosurface(X,Y,Z,reshape(triche(U),size(X)))
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis([0 1 0 1 0 1])
    colorbar
    
    pause
        
    itermov=itermov+1;
    mov(itermov) = getframe(gcf, [0 0 560 420]);
end

%% *** time iterations ****************************************************
while t<=Tmax
    iter=iter+1;
    t=t+ddt;
    clc; disp([iter floor(Tmax/ddt) max(max(U))]);
    
    f=second_membre(X,Y,Z,t,U,U0);
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
    
    b1=-ddt*LW1-ddt*f;
    b2=epsilon*LU1-W+(1/epsilon)*NL;
    [UU,WW] = solver3D(b1, b2, tau, ddt, epsilon);
    U=U+UU;
    W=W+WW;

    if strcmp(film,'yes')==1 && mod(iter,nper)==0
        figure(1)
        isosurface(X,Y,Z,reshape(triche(U),size(X)))
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis([0 1 0 1 0 1])
        
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
isosurface(X,Y,Z,reshape(triche(U),size(X)))
title('corrected image')
xlabel('x')
ylabel('y')
zlabel('z')
grid on
axis([0 1 0 1 0 1])

figure(2)
isosurface(X,Y,Z,reshape(U0,size(X)))
title('perturbated image')
xlabel('x')
ylabel('y')
zlabel('z')
grid on
axis([0 1 0 1 0 1])


figure(3)
isosurface(X,Y,Z,reshape(U,size(X)))
title('image after Cahn-Hilliard equation')
xlabel('x')
ylabel('y')
zlabel('z')
grid on
axis([0 1 0 1 0 1])


figure(4)
isosurface(X,Y,Z,reshape(triche(U)-Unp,size(X)))
title('error')
xlabel('x')
ylabel('y')
zlabel('z')
grid on
axis([0 1 0 1 0 1])
