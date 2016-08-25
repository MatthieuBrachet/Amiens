function U=solv_poiss_neumann_3D(F)
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

  
  %spectre de la matrice
 
  K=alpha+beta*LAPLA;
  fhat=dct(dct(dct(F)')')';
  uhat=fhat./K;
  U=idct(idct(idct(uhat)')')';
 