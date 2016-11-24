function [ SM ] = sm3D( X,Y,Z,t )
global epsilon
u=cos(pi*X).*cos(pi*Y).*cos(pi*Z).*exp(-t);
dtu=-u;
Lu=3*pi^2.*u;
f=u.*(u.^2-1)/epsilon^2;
SM=dtu+Lu+f;
end