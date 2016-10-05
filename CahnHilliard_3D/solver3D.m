function [x1, x2] = solver3D(b1, b2, tau, ddt, epsilon)
%solver for the Cahn-Hilliard system with Schur complement and fft solver.
global A
f=b2+tau*epsilon*A*b1;
alpha=1;
beta=tau^2*epsilon*ddt;
[x2]=prec_bilap_neum3d(f,alpha,beta);
%x2=(Id+tau^2*epsilon*ddt*A^2)\f;
x1=b1-tau*ddt*A*x2;
end