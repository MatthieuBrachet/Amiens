function [f,fnp] = initial_fun(x,y,z)
% 3D initial map.
global test
[n1,n2,n3]=size(x);
if test == 1
    % CI aleatoire pas d'inpainting. CH classique.
    f=2*rand(size(x))-1;
    f=f-mean(mean(mean(f)));
    fnp=f;
    
elseif test == 2
    % bar tachée par une autre bar verticale.
    for i=1:n1
        for j=1:n2
            for k=1:n3
                x1=x(i,j,k); y1=y(i,j,k); z1=z(i,j,k);
                r1=abs(x1-0.6);
                r2=abs(y1-0.6);
                if r1<0.25 & r2<0.15
                    fnp(i,j,k)=1;
                else
                    fnp(i,j,k)=0;
                end
            end
        end
    end
    f=fnp.*(1-(abs(z-0.5)<0.2));
end

