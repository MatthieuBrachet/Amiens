function y=dst(x)
% Discrete Sine Transform DST-I
[M,N]=size(x);
y=[zeros(1,N);x;zeros(1,N);-flipud(x)];
yy=fft(y);
y=real(yy(2:M+1,:)/(-2*i));