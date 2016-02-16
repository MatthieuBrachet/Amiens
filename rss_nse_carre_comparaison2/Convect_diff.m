function [W] = Convect_diff(dt,bx,by,Psi,W,tau,Re)
N=length(Psi); N=sqrt(N);
h=1./(N+1);
[a0,Mx,My,Nx,Ny] = Mlaplacien2(N,2);
%% approximation de -d2W
b=bx+by;
ddW=Mx\(Nx*W+a0*bx)+My\(Ny*W+a0*by);
%% partie non linéaire
[FF] = NL(Psi,W,bx,by);
%% calcul du W suivant
SM=-dt*FF-dt/Re*ddW;
[ZZ] = solvpoiss(SM,N,dt,h,tau,Re);
%ZZ=(Id+dt*tau*A2/Re)\(SM);
W=W+ZZ;
end

