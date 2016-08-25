function [ u ] = sol_exacte(t)
global X Y
u=cos(pi*X).*cos(pi*Y).*exp(-t);
end