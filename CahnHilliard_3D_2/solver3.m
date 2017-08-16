function [ x1,x2 ] = solver3( b1,b2 )
global WW XX YY ZZ
z1=b1-YY*b2;
x1=XX\z1;
x2=b2+ZZ*x1;
end