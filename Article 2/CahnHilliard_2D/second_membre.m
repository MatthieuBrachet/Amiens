function [ f ] = second_membre(x,y,t,U,U0)
global test
if test == 1
    f=0;
    
elseif test == 2
    lambda=100000;
    [n1,n2]=size(x);
    for i=1:n1
        for j=1:n2
            if abs(x(i,j)-0.5)<0.1 & abs(y(i,j)-0.6)<0.3
                xi(i,j)=0;
            else
                xi(i,j)=1;
            end
        end
    end
    xi=reshape(xi,[],1);
    f=lambda.*xi.*(U-U0);
    
elseif test == 3
    lambda=100000;
    [n1,n2]=size(x);
    for i=1:n1
        for j=1:n2
            if abs(x(i,j)-0.5)<0.1 & abs(y(i,j)-0.5)<0.4
                xi(i,j)=0;
            else
                xi(i,j)=1;
            end
        end
    end
    xi=reshape(xi,[],1);
    f=lambda.*xi.*(U-U0);
    
elseif test == 4
    lambda=100000;
    [n1,n2]=size(x);
    for i=1:n1
        for j=1:n2
            if abs(x(i,j)-0.5)<0.25 & abs(y(i,j)-0.5)<0.25
                xi(i,j)=0;
            else
                xi(i,j)=1;
            end
        end
    end
    xi=reshape(xi,[],1);
    f=lambda.*xi.*(U-U0);
    
    
end

end