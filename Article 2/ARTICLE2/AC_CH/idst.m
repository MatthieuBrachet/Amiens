function y=idst(x)
% Inverse Discrete Sine Transform IDST-I
M=size(x,1);
y=2/(M+1)*dst(x);