function [f,fnp] = initial_fun(x,y)
global test
[n1,n2]=size(x);
if test == 1
    f=2*rand(size(x))-1;
    f=f-mean(mean(f));
    fnp=f;
    
elseif test == 2
    fnp=2*(abs(y-0.5)<0.1)-1;
    f=fnp.*(1-(abs(x-0.5)<0.1).*(abs(y-0.6)<0.3));
    
elseif test == 3
    for i=1:n1
        for j=1:n2
            if abs(x(i,j)-0.25)<0.35 & abs(y(i,j)-0.55)<0.1
                fnp(i,j)=1;
            elseif abs(x(i,j)-0.75)<0.35 & abs(y(i,j)-0.4)<0.1
                fnp(i,j)=1;
            else
                fnp(i,j)=-1;
            end
        end
    end    
    f=fnp.*(1-(abs(x-0.5)<0.1).*(abs(y-0.5)<0.4));
    
elseif test == 4
    for i=1:n1
        for j=1:n2
            if (x(i,j)-0.25)^2+(y(i,j)-0.25)^2<0.03
                fnp(i,j)=1;
            elseif (x(i,j)-0.75)^2+(y(i,j)-0.25)^2<0.03
                fnp(i,j)=1;
            elseif (x(i,j)-0.25)^2+(y(i,j)-0.75)^2<0.03
                fnp(i,j)=1;
            elseif (x(i,j)-0.75)^2+(y(i,j)-0.75)^2<0.03
                fnp(i,j)=1;
            else
                fnp(i,j)=-1;
            end
        end
    end    
    f=fnp.*(1-(abs(x-0.5)<0.25).*(abs(y-0.5)<0.25));
    
end
end

