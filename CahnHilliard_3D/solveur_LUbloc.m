function [x1, x2] = solveur_LUbloc(b1, b2, tau, ddt, epsilon)
%solver for the Cahni-Hilliard system with Schur complement
global A Id
x2=(Id+tau^2*epsilon*ddt*A^2)\(b2+tau*epsilon*A*b1);
x1=b1-tau*ddt*A*x2;
end