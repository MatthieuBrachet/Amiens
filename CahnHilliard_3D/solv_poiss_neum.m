function U=solv_poiss_neum(F)
%----------------------------------
%Numerical solution of
%alpha U - \beta \Delta U = F
%
%+ Neumann BC
%
% F must be given as a N^2 column vector 
% u1                 is ---------------
%----------------------------------
global h N alpha beta ; 

  f=reshape(F,N,N);
  %spectre de la matrice
  K=alpha-2*beta*(cos(pi*(0:N-1)'*ones(1,N)*h)+cos(pi*ones(N,1)*(0:N-1)*h)-2)/h^2;
  fhat=dct(dct(f)')';
  uhat=fhat./K;
  u1=idct(idct(uhat)')';
  U=reshape(u1,N*N,1);