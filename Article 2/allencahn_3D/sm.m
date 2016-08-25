function [ SM ] = sm( t )
global X Y Z epsilon
u=cos(pi*X).*cos(pi*Y).*cos(pi*Z).*exp(sin(t));
dtu=cos(t)*u;
Lu=3*pi^2.*u;
SM=dtu+Lu+u.*(u.^2-1)*(1/epsilon^2);
SM=reshape(SM,[],1);
end