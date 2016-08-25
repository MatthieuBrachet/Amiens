function U=solv_poiss_neum_3D(F)
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
global LAPLA ;

  f=reshape(F,N,N,N);
  %spectre de la matrice
 
  K=alpha+beta*LAPLA;
  fhat=dct(dct(dct(f)')')';
  uhat=fhat./K;
  u1=idct(idct(idct(uhat)')')';
  U=reshape(u1,N*N*N,1);