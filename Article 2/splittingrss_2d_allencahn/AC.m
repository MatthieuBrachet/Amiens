function [ Uf, T, E ] = AC(U0, dt, Tmax)
global eps
global A Id
t=0;
E=[]; T=[];
U=U0;
while t<Tmax
    t=t+dt;
    SM=reshape(sm(t),[],1);
    
    F=(1/eps^2).*U.*(U.^2-1);
    Z=(Id+dt*A)\(-dt*A*U-dt*F+dt*SM);
    U=Z+U;

    % calcul de l'erreur
    Uex=reshape(sol_exacte(t),[],1);
    e=norm(U-Uex,inf);
    E=[E e];
    T=[T t];
end
Uf=U;
end

