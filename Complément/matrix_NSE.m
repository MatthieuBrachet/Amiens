%%%%%%%%%%%
%MATRIX_NSE
%%%%%%%%%%%

ALAP=gallery('poisson',N,N)/h^2;
v=[0 1 zeros(1,N-2)];
deix=toeplitz(-v',v)/(2*h);
ident=speye(N,N);
DX=kron(ident,deix);
DY=kron(deix,ident);