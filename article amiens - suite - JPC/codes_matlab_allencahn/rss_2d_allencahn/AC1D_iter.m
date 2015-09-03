function [U] = AC1D_iter(U,f,dt,tau,eps)
% une itération de Allen-Cahn par RSS

    global A2 Mx Nx My Ny Id

    FU=U.*(U.^2-1)./eps^2;
    Ux=Nx*U;
    Uy=Ny*U;
    Z=(Id+dt*tau*A2)\(-dt*(Mx\Ux+My\Uy)+dt*f-dt*FU);
    U=Z+U;
end

