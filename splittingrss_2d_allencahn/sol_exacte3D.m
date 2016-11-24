function [ u ] = sol_exacte3D(X,Y,Z,t)
u=cos(pi*X).*cos(pi*Y).*cos(pi*Z).*exp(-t);
end