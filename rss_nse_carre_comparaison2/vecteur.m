function[F1]=vecteur(F)
    %transforme une matrice en vecteur
    F=F';
    F1=[];
    for i=1:size(F,1)
        F1=[F1;F(:,i)];
    end
end