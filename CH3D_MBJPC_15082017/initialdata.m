function [ u0, u, ind ] = initialdata(X,Y,Z,test)

if strcmp(test,'box')==1
    ind=1-(abs(X-.5)<.45).*(abs(Y-.45)<.45).*(abs(Z-.5)<.07);
    u0=2*(abs(X-.5)<.1).*(abs(Y-.5)<.1).*(abs(Z-.5)<.49)-1;
    u0=u0.*ind;
    u=u0;
elseif strcmp(test,'scroll')==1
    ind=1-(abs(X-.75)<.22).*(abs(Y-.45)<.025).*(abs(Z-.5)<.1);
    
    Gradius=0.3;
    Pradius=0.2;
    a=4*pi;

    [n1,n2,n3]=size(X);
    for i=1:n1
        for j=1:n2
            for k=1:n3
                xx=X(i,j,k);
                yy=Y(i,j,k);
                zz=Z(i,j,k);
                M=[Gradius.*cos(a*zz)+.5; Gradius.*sin(a*zz)+.5; zz];
                P=[xx;yy;zz];
                dist=norm(M-P,2);
                f(i,j,k)=(dist<Pradius);
            end
        end
    end
    u0=f.*ind;
    u=u0;
end

end

