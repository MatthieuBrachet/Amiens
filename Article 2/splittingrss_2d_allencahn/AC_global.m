function [ Uf, T, E ] = AC_global(U0, dt, Tmax)
global eps
global A Id 
t=0;
E=[]; T=[];
U=U0;
while t<Tmax
    t=t+dt;
    
    % resolution de la partie non linéaire
    ddt=dt/2;
    nom=U;
    denom=U.^2+(1-U.^2).*exp(-2*ddt/(eps^2));
    denom=sqrt(denom);
    U=nom./denom;
    
    % resolution de l'equation de la chaleur
    ddt=dt;
    SM=reshape(sm(t),[],1);
    Z=(Id+ddt*A)\(-ddt*A*U + ddt*SM);
    U=Z+U;
    
    % resolution de la partie non linéaire
    ddt=dt/2;
    nom=U;
    denom=U.^2+(1-U.^2).*exp(-2*ddt/(eps^2));
    denom=sqrt(denom);
    U=nom./denom;

    % calcul de l'erreur
    Uex=reshape(sol_exacte(t),[],1);
    e=norm(U-Uex,inf)./norm(Uex,inf);
    E=[E e];
    T=[T t];
end
Uf=U;
end