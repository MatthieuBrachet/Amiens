function [ u ] = sol_exacte(t)
global X Y Z
u=cos(pi*X).*cos(pi*Y).*cos(pi*Z).*exp(sin(t));
end