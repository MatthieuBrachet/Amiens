function y = dst(x) 
n = size(x,1); m = size(x,2);
y = [zeros(1,m);x]; 
y = imag(fft(y,2*n+2)); 
y = sqrt(2/(n+1))*y(2:n+1,:);
return