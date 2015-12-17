function [AA] = level(A,a,b)
%fonction pour extraire les valeurs de A comprises entre a et b

AA=zeros(size(A));

%% supression de la partie supérieure et inférieure
for i=1:size(A,1)
    for j=1:size(A,2)
        if a <= A(i,j) & A(i,j) <= b
            AA(i,j)=A(i,j);
        else
            AA(i,j)=0;
        end
    end
end

end

