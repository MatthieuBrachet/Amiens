function [ut] = triche(u)
% remplace les valeurs de u par 1 ou -1 en fonction de >0 ou <0 (0 est maintenu). 
[n1,n2]=size(u);
for i=1:n1
    for j=1:n2
        if u(i,j)>0
            ut(i,j)=1;
        elseif u(i,j)<0
            ut(i,j)=-1;
        else
            ut(i,j)=0;
        end
    end
end
end