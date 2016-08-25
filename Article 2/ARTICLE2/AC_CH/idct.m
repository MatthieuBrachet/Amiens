function y=idct(x)
% Inverse Discrete Cosine Transform IDCT-I
M=size(x,1);
y=2/(M-1)*dct(x);