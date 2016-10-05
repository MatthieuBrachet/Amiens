function y=dct(x)
% Discrete Cosine Transform DCT-I
[M,N]=size(x);
y=[x;flipud(x(2:M-1,:))];
yy=fft(y);
y=real(yy(1:M,:)/2);