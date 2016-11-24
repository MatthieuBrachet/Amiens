function [y] = indicatrice(x,y,a,b)
ind=(abs(x-.5)<a/2).*(abs(y-.5)<b/2);
y=1-ind;
end

