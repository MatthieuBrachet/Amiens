function [W] = Convect_diff(dt,bx,by,Psi,W,tau)
global N Re h Mx My Nx Ny a0
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

