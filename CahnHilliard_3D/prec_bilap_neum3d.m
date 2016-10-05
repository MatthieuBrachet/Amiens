function [u]=prec_bilap_neum3d(f,alpha,beta)
%% *****************************************
% solve :
%
%    alpha * u + \Delta^2 u = f
%
% f must be a N*N*N vector.
% u is a N*N*N vector.

global nb h
one=ones(nb,1);
K=cos(pi*(0:nb-1)*h);
%
DX2=kron(one,K);
DY2=DX2';


for j=1:nb
     DDX2(:,:,j)=2*(1-DX2);
     DDY2(:,:,j)=2*(1-DY2);
     DDZ2(j,:,:)=2*(1-DX2);
end
f=reshape(f,nb,nb,nb);  
LAPLA=(DDX2+DDY2+DDZ2)/h^2;
K=alpha+beta*LAPLA.^2;
fc=mirt_dctn(mirt_dctn(f));
uc=fc./K;
U=mirt_idctn(mirt_idctn(uc));
u=reshape(U,[],1);