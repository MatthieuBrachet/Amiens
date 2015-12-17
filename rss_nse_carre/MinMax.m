function [MIN,MAX] = MinMax(PP,X,Y)
%SORTIE :
%     MIN = (Min(PP), x_min, y_min)
%     MAX = (Max(PP), x_max, y_max)
%
%ENTREE : 
%     PP : courbe
%     X  : coordonn�es en x. X est une matrice de m�me taille que PP.
%     Y  : coordonn�es en y. Y est une matrice de m�me taille que PP.
%     X et Y peuvent �tre construits avec meshgrid.


%% MAX

Max=max(max(PP));
for i=1:size(PP,1)
    for j=1:size(PP,2)
        if PP(i,j) == Max
            MAX=[Max, X(i,j), Y(i,j)];
        end
    end
end

%% MIN

Min=min(min(PP));
for i=1:size(PP,1)
    for j=1:size(PP,2)
        if PP(i,j) == Min
            MIN=[Min, X(i,j), Y(i,j)];
        end
    end
end


end