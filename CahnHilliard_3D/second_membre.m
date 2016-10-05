function [ f ] = second_membre(x,y,z,t,U,U0)
% 3D r.h.s. map
global test
[n1,n2,n3]=size(x);
lambda=500000;
if test == 1
    f=0;
elseif test == 2
    for i=1:n1
        for j=1:n2
            for k=1:n3
                if abs(z(i,j,k)-0.5)<0.2
                    xi(i,j,k)=0;
                else
                    xi(i,j,k)=1;
                end
            end
        end
    end
    xi=reshape(xi,[],1);
    f=lambda.*xi.*(U-U0);
end

end