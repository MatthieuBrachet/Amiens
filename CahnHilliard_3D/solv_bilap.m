function U=solv_bilap(F)
%--------------------------------------
%Numerical solution of
%alpha U + \beta \Delta^2 U = F
%
%+ Neumann BC
%
% F must be given as a N^3 column vector 
% u1 is --------------------------------

global h N alpha beta K

DX2=cos(pi*(0:N-1)'*ones(1,N)*h);
DY2=cos(pi*ones(N,1)*(0:N-1)*h);
for j=1:N
    DDX2(:,:,j)=2*(1-DX2);
    DDY2(:,:,j)=2*(1-DY2);
end
for k=1:N
    DDZ2(k,:,:)=2*(1-DX2);
end
      
LAPLA=(DDX2+DDY2+DDZ2)/h^2;

K=alpha + beta * LAPLA.^2;
Fhathat=dctn(dctn(F));
Uhathat=Fhathat./K;
U=idctn(idctn(Uhathat));


