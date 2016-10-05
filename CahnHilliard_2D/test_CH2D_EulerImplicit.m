clc; clear all; close all;

%% ************************************************************************
global test

%% *** options ************************************************************
film='yes';
test=1;

%% *** numerical data
n=40;
[a0,id,AA] = Mlaplacien2(n,2);
Ax=kron(AA,id);
Ay=kron(id,AA);
A=Ax+Ay;
Id=speye(size(A));

%% *** space data *********************************************************
h=1/(n+1);
x=[0:h:1]';
[X,Y]=meshgrid(x,x);

%% *** time data **********************************************************
t=0;
ddt=10^-6;
Tmax=0.1;

%% *** mathematical data **************************************************
epsilon=0.01;
u=initial_fun(X,Y);
U=reshape(u,[],1);
W=epsilon*A*U+(1/epsilon).*U.*(U.^2-1);

MAT=[Id ddt*A; -epsilon*A Id];
%% *** time iterations ****************************************************
while t<Tmax
    clc; t=t+ddt
    
    b=[-ddt*A*W; epsilon*A*U-W+(1/epsilon)*U.*(U.^2-1)];
    Z=MAT\b;
    
    U=U+Z(1:(n+2)*(n+2),1);
    W=W+Z((n+2)*(n+2)+1:end,1);
    
    if strcmp(film,'yes')==1
        % pause
    
        figure(1)
        contourf(X,Y,reshape(U,size(X)))
        xlabel('x')
        ylabel('y')
        colorbar
    end
end
   
figure(1)
contourf(X,Y,reshape(U,size(X)))
xlabel('x')
ylabel('y')
colorbar





