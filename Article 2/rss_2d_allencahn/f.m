function [F] = f(X,Y,t,eps)
%second membre pour AC
U=cos(pi*X).*cos(pi*Y)*exp(sin(t));
F=(cos(t)+2*pi^2).*U+U.*(U.^2-1)*1/(eps^2);
F=reshape(F,[],1);
end

