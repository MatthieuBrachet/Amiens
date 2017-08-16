function [ u0, u, ind ] = initialdata(X,Y,test)

if strcmp(test,'box')==1
    ind=1-(abs(X-.5)<.3).*(abs(Y-.5)<.05);
    u0=2*(abs(X-.5)<.1)-1;
    u0=u0.*ind;
    u=u0;
    
elseif strcmp(test,'circles')==1
    ind=1-(abs(X-.5)<.25).*(abs(Y-.5)<.25);
    rayon=0.1;
    u0=(sqrt((X-.25).^2+(Y-.25).^2)<rayon);
    u0=u0+(sqrt((X-.75).^2+(Y-.25).^2)<rayon);
    u0=u0+(sqrt((X-.25).^2+(Y-.75).^2)<rayon);
    u0=u0+(sqrt((X-.75).^2+(Y-.75).^2)<rayon);
    u0=2*u0-1;
    u0=u0.*ind;
    u=u0;
    
elseif strcmp(test,'triangle')==1
    ind=1-(X>=0.02).*(X<=0.98).*(Y>=0.34).*(Y<=.44);
    f1= (Y>=0.2).*(Y<=2*X-.2).*(Y<=-2*X+1.8);
    f1=2*f1-1;
    u0=f1.*ind;
    u=u0;
end

end

