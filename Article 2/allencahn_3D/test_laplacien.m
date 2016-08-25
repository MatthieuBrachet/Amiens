clc; clear all; close all;

NN=10:10:100

for i=1:length(NN)

    clc; N=NN(i)
    h=1/(N+1);
    x=[0:h:1]';
    [X,Y,Z]=meshgrid(x,x,x);

    A=sparse(diag(-2*ones(N+2,1))+diag(ones(N+1,1),1)+diag(ones(N+1,1),-1));
    A(1,2)=2; A(end,end-1)=2;
    A=-A./(h^2);
    id=speye(N+2,N+2);

    Ax=kron(kron(A,id), id);
    Ay=kron(kron(id,A), id);
    Az=kron(kron(id,id), A);

    A=Ax+Ay+Az;

    %% mathématical data
    u=cos(pi*X).*cos(pi*Y).*cos(pi*Z);
    u=reshape(u,[],1);

    f=3*pi^2*u;

    ee(i)=max(abs(A*u-f));
    hh(i)=h;
end

figure(1)
loglog(hh,ee,hh,hh.^2)
    