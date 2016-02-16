function [Psi,W0] = Stokes(N,Re)
h=1./(N+1);
x=[h:h:1-h];

%matrice
A=sparse(gallery('poisson',N)/h^2);
W0=kron(sin(pi*x),sin(pi*x))';

kmax=20;
k=0;
res=10;
RESIDU=[];
while k< kmax & res >1.e-8;
    Psi=-A\(W0);
    [bwx,bwy] = bord(Psi);
    F=(bwx+bwy);
    W0=Re*A\F;
    k=k+1;
end
end

