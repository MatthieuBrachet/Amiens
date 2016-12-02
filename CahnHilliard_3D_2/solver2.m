function [ x1,x2 ] = solver2( b1,b2 )
global WW XX YY ZZ
z1=b1;
z2=b2-XX*z1;
x2=WW\z2;
x1=YY\(z1-ZZ*x2);
end