function [U]=prec_bilap_neum3d(f,alpha,beta)
global N h
one=ones(N,1);
K=cos(pi*(0:N-1)*h);
%
DX2=kron(one,K);
DY2=DX2';


for j=1:N
     DDX2(:,:,j)=2*(1-DX2);
     DDY2(:,:,j)=2*(1-DY2);
     DDZ2(j,:,:)=2*(1-DX2);
end
f=reshape(f,N,N,N);  
LAPLA=(DDX2+DDY2+DDZ2)/h^2;
K=alpha+beta*LAPLA.^2;
fc=mirt_dctn(f);
uc=fc./K;
U=mirt_idctn(uc);