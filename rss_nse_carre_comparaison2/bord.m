function [bwx,bwy] = bord(psi)
%mise à jour des conditions de bord
global cavite

N=length(psi); N=sqrt(N);
h=1/(N+1);
x=h:h:1-h;

A=speye(N,N);
B=zeros(N,N);
B(1,1)=8;
B(1,2)=-3;
B(1,3)=8/9;
B(1,4)=-1/8;

B(end,end)=B(1,1);
B(end,end-1)=B(1,2);
B(end,end-2)=B(1,3);
B(end,end-3)=B(1,4);

Kx=(1/(h^2))*kron(A,B);
Ky=(1/(h^2))*kron(B,A);

b=zeros(N,N);
if cavite == 1
    b(end,:)=1;
elseif cavite == 2
    b(end,:)=(1-(1-2*x).^2).^2;
else
    error('error cavity definition.')
end
bb=vecteur(b);

bwx=Kx*psi;
bwy=Ky*psi+25/(6*h)*bb;
end

