function [a0,Mx,My,Nx,Ny] = Mlaplacien2(n,dimension)
%approximation compacte de '-laplacien' a l'ordre 4
%
% N : nombre de discr�tisation en espace
% dimension : dimension de l'espace de travail (1 ou 2)
%
    h=1/(n+1);
    
    Id=speye(n,n);
    
    M=speye(n,n)+sparse(1/10*diag(ones(n-1,1),1)+1/10*diag(ones(n-1,1),-1));
    
    N=2*sparse(eye(n,n))-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);
    N=6/5*N;
    N(1,1)=67/60;
    N(1,2)=7/12;
    N(1,3)=-13/10;
    N(1,4)=61/120;
    N(1,5)=-1/12;
    N(end,end)=N(1,1);
    N(end,end-1)=N(1,2);
    N(end,end-2)=N(1,3);
    N(end,end-3)=N(1,4);
    N(end,end-4)=N(1,5);
    N=sparse((1/h^2)*N);
    
    a0=-33/40*(1/h^2);
    
    if dimension == 1
        Mx=M;
        Ny=N;
        a0=a0;
        My=M;
        Nx=N;
        
    elseif dimension == 2
        
        Mx=kron(Id,M);
        My=kron(M,Id);
        Nx=kron(Id,N);
        Ny=kron(N,Id);
        a0=a0;
end

