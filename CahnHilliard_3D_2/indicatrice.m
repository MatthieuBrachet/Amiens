function [y] = indicatrice(x,y,z,a,b,c)
global centre
ind=(abs(x-centre(1))<a/2).*(abs(y-centre(2))<b/2).*(abs(z-centre(3))<c/2);
y=1-ind;
end

