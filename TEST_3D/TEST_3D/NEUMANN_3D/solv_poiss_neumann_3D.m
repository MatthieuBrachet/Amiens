function U=solv_poiss_neumann_3D(F,alpha,beta)
%----------------------------------
%Numerical solution of
%alpha U - \beta \Delta U = F
%
%+ Neumann BC
%
% F must be given as a N^2 column vector 
% u1                 is ---------------
%----------------------------------
[N,Nwk1,Nwk2]=size(F);
h=1/(N-1);
%% 3D Neumann Frequency Matrix
one=ones(N,1);
K=cos(pi*(0:N-1)*h);
%
DX2=kron(one,K);
DY2=DX2';

%DX2=cos(pi*(0:N-1)'*ones(1,N)*h);
%DY2=cos(pi*ones(N,1)*(0:N-1)*h);
for j=1:N
     DDX2(:,:,j)=2*(1-DX2);
     DDY2(:,:,j)=2*(1-DY2);
     DDZ2(j,:,:)=2*(1-DX2);
end
      
LAPLA=(DDX2+DDY2+DDZ2)/h^2;
K=alpha+beta*LAPLA;

%% solveur
UC=mirt_dctn(F)./K;
U=mirt_idctn(UC);