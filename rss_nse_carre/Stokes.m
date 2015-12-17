function [Psi,W0] = Stokes()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global x N h Re
%matrice
A=sparse(gallery('poisson',N)/h^2);
W0=kron(sin(pi*x),sin(pi*x));

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

