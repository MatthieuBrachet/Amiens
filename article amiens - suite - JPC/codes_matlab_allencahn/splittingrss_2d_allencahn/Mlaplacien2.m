function [a0,Mx,Nx] = Mlaplacien2(n)
%approximation compacte de '-laplacien' a l'ordre 4 + neumann homogène
%
% n : nombre de discrï¿½tisation en espace
%
    h=1/(n+1);

    %% axe des x
    Mx=speye(n+2,n+2)+sparse(1/10*diag(ones(n+1,1),1)+1/10*diag(ones(n+1,1),-1));
    
    N=2*sparse(eye(n+2,n+2))-diag(ones(n+1,1),1)-diag(ones(n+1,1),-1);
    N=6/5*N;
    N(1,1)=2681/480;
    N(1,2)=-23/3;
    N(1,3)=113/40;
    N(1,4)=-13/15;
    N(1,5)=59/480;
    N(end,end)=N(1,1);
    N(end,end-1)=N(1,2);
    N(end,end-2)=N(1,3);
    N(end,end-3)=N(1,4);
    N(end,end-4)=N(1,5);
    Nx=sparse(N)/(h^2);

    %% pour le second membre
    a0=-33/40*(1/h^2);
    
end
