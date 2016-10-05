clc; clear all; close all;

%% ************************************************************************
global test

%% *** options ************************************************************
film='yes';
test=1;

%% *** numerical data
n=20;
[a0,id,AA] = Mlaplacien2(n,2);
Ax=kron(kron(AA,id),id);
Ay=kron(kron(id,AA),id);
Az=kron(kron(id,id),AA);
A=Ax+Ay+Az;
Id=speye(size(A));

%% *** space data *********************************************************
h=1/(n+1);
x=[0:h:1]';
[X,Y,Z]=meshgrid(x,x,x);

%% *** time data **********************************************************
t=0;
ddt=10^-6;
Tmax=0.01;

%% *** mathematical data **************************************************
epsilon=0.01;
u=initial_fun(X,Y,Z);
U=reshape(u,[],1);
W=epsilon*A*U+(1/epsilon).*U.*(U.^2-1);

MAT=[Id ddt*A; -epsilon*A Id];
%% *** time iterations ****************************************************
while t<Tmax
    clc; t=t+ddt
    
    b=[-ddt*A*W; epsilon*A*U-W+(1/epsilon)*U.*(U.^2-1)];
    ZZ=MAT\b;
    
    U=U+ZZ(1:(n+2)*(n+2)*(n+2),1);
    W=W+ZZ((n+2)*(n+2)*(n+2)+1:end,1);
    
    if strcmp(film,'yes')==1
        % pause
    
        figure(1)
        isosurface(X,Y,Z,reshape(U,size(X)))
        xlabel('x')
        ylabel('y')
        zlabel('z')
        hold off
        close(1)
    end
end



figure(1)
isosurface(X,Y,Z,reshape(U,size(X)))
xlabel('x')
ylabel('y')
zlabel('z')





