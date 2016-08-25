clc; clear all; close all;

NN=10:10:1000;
order=4;
for i=1:length(NN)
    N=NN(i);
    h=1/(N+1);  
    x=[0:h:1]';

    u=cos(pi*x);

    [a0,M,N] = Mlaplacien2(N,order);
    
    

    ddu=(M\N)*u;    
    dduex=pi^2*u;

    E(i)=max(abs(ddu-dduex));
    H(i)=h;
end

figure(1)
semilogy(x,abs(ddu-dduex))

figure(2)
plot(x,ddu,x,dduex,'r-')

figure(3)
loglog(H,E,H,H.^order)