function [f,fp] = initial_fun(x,y,z)
global barx bary barz
global test

if test ==1
    f=(abs(x-0.5)<0.1).*(abs(y-0.5)<1).*(abs(z-0.5)<.1);
elseif test == 2
    Gradius=0.3;
    Pradius=0.2;
    a=4*pi;

    [n1,n2,n3]=size(x);
    for i=1:n1
        for j=1:n2
            for k=1:n3
                xx=x(i,j,k);
                yy=y(i,j,k);
                zz=z(i,j,k);
                M=[Gradius.*cos(a*zz)+.5; Gradius.*sin(a*zz)+.5; zz];
                P=[xx;yy;zz];
                dist=norm(M-P,2);
                f(i,j,k)=(dist<Pradius);
            end
        end
    end
elseif test == 3
    f=(abs(x-0.5)<0.1).*(abs(y-0.5)<1).*(abs(z-0.5)<.1);
    f=f+(abs(x-0.25)<0.1).*(abs(y-0.25)<1).*(abs(z-0.25)<.1);
    f=f+(abs(x-0.75)<0.1).*(abs(y-0.75)<1).*(abs(z-0.75)<.1);
end
f=2*f-1;
ind=indicatrice(x,y,z,barx,bary,barz);
fp=f.*(ind);