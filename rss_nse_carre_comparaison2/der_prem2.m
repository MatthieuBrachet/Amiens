function [a0,Mx,My,Nx,Ny] = der_prem2(n,dimension)
%
h=1/(n+1);
N=diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
N=3/2*N;
N(1,1)=-2;
N(1,2)=3;
N(1,3)=-2/3;
N(1,4)=1/8;
            
N(end,end)=2;
N(end,end-1)=-3;
N(end,end-2)=2/3;
N(end,end-3)=-1/8;
            
M=1/4*diag(ones(n-1,1),1)+1/4*diag(ones(n-1,1),-1)+speye(n,n);
            
N=(1/(2*h))*N;
         
a0=-11/24;  

Id=speye(n);
if dimension == 1
    Mx=M;
    Ny=N;
    a0=a0/(2*h);
    My=M;
    Nx=N;
    
elseif dimension == 2
    
    Mx=sparse(kron(Id,M));
    My=sparse(kron(M,Id));
    Nx=sparse(kron(Id,N));
    Ny=sparse(kron(N,Id));
    a0=a0/(2*h);

end

