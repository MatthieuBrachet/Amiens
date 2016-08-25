function [gh,ghp,lambh]=initial(s,t)
global initial_type ;
lamb= 90000;
switch initial_type
    case 'cercles'
[x1,y1]=meshgrid(s,t);
x=x1/2;y=y1/2;
f1=  ((x-0.12).^2 + (y-0.12).^2 <= 0.08^2)...
+((x-0.37).^2 + (y-0.37).^2 <= 0.08^2) + ...
((x-0.12).^2 + (y-0.37).^2 <= 0.08^2)...
+((x-0.37).^2 + (y-0.12).^2 <= 0.08^2);

f2=(y>=0.12).*(y<=0.37).*(x>=0.12).*(x<=0.37);


%func f3=f1*(1-f2);
%func f4=((x-0.25)^2+(y-0.25)^2 > 0.15^2)*(1 -f2);
%gh=f3-f2-f4;
gh=f1.*(1-f2);

ghp=gh+ rand(1)*(y>=0.12).*(y<=0.37).*(x>=0.12).*(x<=0.37);

lambh=lamb*(1-(y>=0.12).*(y<=0.37).*(x>=0.12).*(x<=0.37))+ 0.*(y>=0.12).*(y<=0.37).*(x>=0.12).*(x<=0.37);
    case 'triangles'
      [x1,y1]=meshgrid(s,t);
      x=x1/2;y=y1/2;
        f1= (y>=0.1).*( y <= 2*x - 0.1).*(y<= -2*x+0.9);
        f2=(y>=0.17).*(y<=0.22);
        
        gh=f1.*(1-f2);


        ghp=gh+ rand(1)*(x>=0.01).*(x<=0.49).*(y>=0.17).*(y<=0.22);
        lambh=lamb*(1-(x>=0.01).*(x<=0.49).*(y>=0.17).*(y<=0.22))+ 0.*(x>=0.01).*(x<=0.49).*(y>=0.17).*(y<=0.22);
end