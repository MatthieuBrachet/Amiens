function [W] = Convect_diff(dt,bx,by,Psi,W,tau,Re)
global RSS


N=sqrt(length(Psi));
h=1/(N+1);
[a0,Mx,My,Nx,Ny] = Mlaplacien2(N,2);
Id=speye(N,N);
%% approximation de -d2W
b=bx+by;
ddW=Mx\(Nx*W+a0*bx)+My\(Ny*W+a0*by);




%% partie non linéaire
[FF] = NL(Psi,W,bx,by);


%% calcul du W suivant
SM=-dt*FF-dt/Re*ddW;
switch  RSS
    case'lineaire'
        [ZZ] = solvpoiss(SM,N,dt,h,tau,Re);
    case 'nonlineaire'
        ALAP=gallery('poisson',N,N)/h^2;
        v=[0 1 zeros(1,N-2)];
        deix=toeplitz(-v',v)/(2*h);
        ident=speye(N,N);
        DX=kron(ident,deix);
        DY=kron(deix,ident);
        ID=speye(size(ALAP));
        
        MATRIX=ID+tau*dt*(ALAP/Re+(sparse(diag(DY*Psi))*DX-sparse(diag(DX*Psi))*DY));
        ZZ=MATRIX\SM;
end
%ZZ=(Id+dt*tau*A2/Re)\(SM);
W=W+ZZ;

end

