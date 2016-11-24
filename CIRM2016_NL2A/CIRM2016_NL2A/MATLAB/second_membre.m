function [ f ] = second_membre(x,y,t,U,U0)
global test

lambda=90000;
[n1,n2]=size(x);
if test == 1
    f=0;
    
elseif test == 2
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
    
elseif test == 5
    for i=1:n1
        for j=1:n2
            if abs(x(i,j)-0.5)<0.1 & abs(y(i,j)-0.5)<0.45
                xi(i,j)=0;
            else
                xi(i,j)=1;
            end
        end
    end
    xi=reshape(xi,[],1);
    f=lambda.*xi.*(U-U0);
    
elseif test == 6
    rec1=[0.25 0.5];
    rec2=[0.5 0.5];
    rec3=[0.75 0.5];
    for i=1:n1
        for j=1:n2
            if abs(x(i,j)-rec1(1))<0.015 & abs(y(i,j)-rec1(2))<0.45
                xi(i,j)=0;
            elseif abs(x(i,j)-rec2(1))<0.015 & abs(y(i,j)-rec2(2))<0.45
                xi(i,j)=0;
            elseif abs(x(i,j)-rec3(1))<0.015 & abs(y(i,j)-rec3(2))<0.45
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