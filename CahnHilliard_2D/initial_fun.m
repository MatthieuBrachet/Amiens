function [f,fnp] = initial_fun(x,y)
global test
[n1,n2]=size(x);
if test == 1
    % CI aleatoire pas d'inpainting. CH classique.
    f=2*rand(size(x))-1;
    f=f-mean(mean(f));
    fnp=f;
    
elseif test == 2
    % bar tachée par une autre bar verticale.
    fnp=2*(abs(y-0.5)<0.1)-1;
    f=fnp.*(1-(abs(x-0.5)<0.1).*(abs(y-0.6)<0.3));
    
elseif test == 3
    % deux bar non alignée tachée d'une bar verticale.
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
    % 4 cercles tachées d'un carré.
    radius=0.15;
    for i=1:n1
        for j=1:n2
            if (x(i,j)-0.25)^2+(y(i,j)-0.25)^2<radius^2
                fnp(i,j)=1;
            elseif (x(i,j)-0.75)^2+(y(i,j)-0.25)^2<radius^2
                fnp(i,j)=1;
            elseif (x(i,j)-0.25)^2+(y(i,j)-0.75)^2<radius^2
                fnp(i,j)=1;
            elseif (x(i,j)-0.75)^2+(y(i,j)-0.75)^2<radius^2
                fnp(i,j)=1;
            else
                fnp(i,j)=-1;
            end
        end
    end    
    f=fnp.*(1-(abs(x-0.5)<0.25).*(abs(y-0.5)<0.25));
    
elseif test == 5
    % un cercle taché d'une bar.
     for i=1:n1
        for j=1:n2
            if (x(i,j)-0.5)^2+(y(i,j)-0.5)^2<0.35^2
                fnp(i,j)=1;
            else
                fnp(i,j)=-1;
            end
        end
    end    
    f=fnp.*(1-(abs(x-0.5)<0.1).*(abs(y-0.5)<0.45));
    
elseif test == 6
    % cosinus a corriger
    for i=1:n1
        for j=1:n2
            if abs(y(i,j)-(0.3*sin(4*pi*x(i,j))+0.5))<0.06
                fnp(i,j)=1;
            else
                fnp(i,j)=-1;
            end
        end
    end
    rec1=[0.25 0.5];
    rec2=[0.5 0.5];
    rec3=[0.75 0.5];
    f=fnp.*(1-(abs(x-rec1(1))<0.015).*(abs(y-rec1(2))<0.45)).*(1-(abs(x-rec2(1))<0.015).*(abs(y-rec2(2))<0.45)).*(1-(abs(x-rec3(1))<0.015).*(abs(y-rec3(2))<0.45));
                
end
end

