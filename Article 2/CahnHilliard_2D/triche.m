function [ut] = triche(u)
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