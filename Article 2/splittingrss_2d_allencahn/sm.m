function [ SM ] = sm( t )
global X Y eps
u=cos(pi*X).*cos(pi*Y).*exp(-t);
dtu=-u;
Lu=2*pi^2.*u;
f=u.*(u.^2-1)/eps^2;
SM=dtu+Lu+f;
end