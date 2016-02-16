function[U1]=matrice(U)
    %transforme un vecteur en matrice
    n=sqrt(length(U));
    for i=1:n
        U1([1:n],i)=U([(i-1)*n+1:i*n]);
    end
    U1=U1';
end