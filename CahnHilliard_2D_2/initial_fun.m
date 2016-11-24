function [f,fp] = initial_fun(x,y)
global barx bary
global test

if test ==1
    f=(abs(x-0.5)<0.1).*(abs(y-0.5)<1);
elseif test ==2
    r1=sqrt((x-.25).^2+(y-.25).^2);
    r2=sqrt((x-.25).^2+(y-.75).^2);
    r3=sqrt((x-.75).^2+(y-.75).^2);
    r4=sqrt((x-.75).^2+(y-.25).^2);
    radius=.2;
    f=(r1<radius)+(r2<radius)+(r3<radius)+(r4<radius);
end
f=2*f-1;
ind=indicatrice(x,y,barx,bary);
fp=f.*(ind);