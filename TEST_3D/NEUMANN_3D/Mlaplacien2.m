function [a0,M,N] = Mlaplacien2(n,order)
%approximation compacte de '-Delta' a l'ordre 4 + neumann homog�ne
% n : nombre de discr�tisation en espace
% M : implicit part,
% N : explicit part.
    h=1/(n+1);

if order == 4
    M=eye(n+2,n+2)+1/10*diag(ones(n+1,1),1)+1/10*diag(ones(n+1,1),-1);
    M=sparse(M);
    
    N=(-12/5*eye(n+2,n+2)+6/5*diag(ones(n+1,1),1)+6/5*diag(ones(n+1,1),-1));
    N(1,1)=-2681/480;
    N(1,2)=23/3;
    N(1,3)=-113/40;
    N(1,4)=13/15;
    N(1,5)=-59/480;
    N(end,end)=N(1,1);
    N(end,end-1)=N(1,2);
    N(end,end-2)=N(1,3);
    N(end,end-3)=N(1,4);
    N(end,end-4)=N(1,5);
    N=-N/(h^2);
    a0=33/(40*h*h);
else
    a0=1;
     M=speye(n+2,n+2);
% %     N=2*sparse(eye(n+2,n+2))-diag(ones(n+1,1),1)-diag(ones(n+1,1),-1);
% %     N(1,1)=1; N(end,end)=1;
% %     N=N./(h^2);
%     a0=1;
%     
%     
%      A=sparse(diag(-2*ones(n+2,1))+diag(ones(n+1,1),1)+diag(ones(n+1,1),-1));
%      A(1,2)=2; A(end,end-1)=2;
%      N=-A./(h^2);
     DX=gallery('tridiag',n+2);
     DX(1,1)=1;DX(n+2,n+2)=1;
     N=DX/h^2;
end