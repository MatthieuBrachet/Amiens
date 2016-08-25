function[a0,A2]=Mlaplacien(N)
    % n : discretisation i.e. dimension de la matrice
    % dicrÃ©tisation laplacien 2d a l'ordre 2
    % neumann homogène

    h=1/(N+1);
    A2=2*speye(N+2,N+2)-1*diag(ones(N+1,1),1)-1*diag(ones(N+1,1),-1);
    A2(1,1)=1;
    A2(1,2)=-1;

    A2(end,end)=A2(1,1);
    A2(end,end-1)=A2(1,2);

    A2=A2./h^2;
    a0=1/h^2;
end