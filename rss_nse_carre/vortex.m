function [ PP1, PP2, PP3, PP4, PP5, PP6, PP7, PP8, PP9 ] = vortex(PP)
    [n,m]=size(PP);
    
    PP1=PP(1:floor(n/3),1:floor(n/3));
    PP2=PP(1:floor(n/3),floor(n/3)+1:2*floor(n/3));
    PP3=PP(1:floor(n/3),2*floor(n/3)+1:end);
    
    PP4=PP(floor(n/3)+1:2*floor(n/3),1:floor(n/3));
    PP5=PP(floor(n/3)+1:2*floor(n/3),floor(n/3)+1:2*floor(n/3));
    PP6=PP(floor(n/3)+1:2*floor(n/3),2*floor(n/3)+1:end);
    
    PP7=PP(2*floor(n/3)+1:end,1:floor(n/3));
    PP8=PP(2*floor(n/3)+1:end,floor(n/3)+1:2*floor(n/3));
    PP9=PP(2*floor(n/3)+1:end,2*floor(n/3)+1:end);
end

