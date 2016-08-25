function [ U ] = nonlinear_part( U, dt )
% non linear part of Allen-Cahn equation solver
global epsilon;

nom=U;
denom=U.^2+(1-U.^2).*exp(-2*dt/epsilon^2);
denom=sqrt(denom);

U=nom./denom;
end

