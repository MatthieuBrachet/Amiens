function [ut] = triche(u)
% remplace les valeurs de u par 1 ou -1 en fonction de >0 ou <0 (0 est maintenu). 
[n1,n2,n3]=size(u);
for i=1:n1
    for j=1:n2
        for k=1:n3
            if u(i,j,k)>0.5
                ut(i,j,k)=1;
            elseif u(i,j,k)<0.5
                ut(i,j,k)=0;
            else
                ut(i,j,k)=0;
            end
        end
    end
end
end