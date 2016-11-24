function [x1, x2] = solveur_LUbloc(b1, b2, tau, ddt, epsilon)
%solver for the Cahn-Hilliard system with Schur complement
global A Id KK
global NX
global label
switch label
    case 'fft'
bp=reshape(b2+tau*epsilon*A*b1,NX,NX);

scmb=dctn(bp);
xc2=idctn(scmb./KK);
x2=reshape(xc2,NX*NX,1);
x1=b1-tau*ddt*A*x2;
    case 'classic'
x2=(Id+tau^2*epsilon*ddt*A^2)\(b2+tau*epsilon*A*b1);
x1=b1-tau*ddt*A*x2;
end


end