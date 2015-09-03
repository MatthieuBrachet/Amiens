function [ U ] = chaleur( U , dt)
    global Ax Ay I Mx Nx My Ny taux tauy
    
    % solver for heat equation (Neumann boundaries condition)
    % with strang splitting
    
    % dt/2 in x-axes
    Ux=Nx*U;
    W=Mx\Ux;
    Z=(I+taux*dt/2*Ax)\(-dt/2*W);
    U1=Z+U;
    
    % dt in y-axes
    Uy=Ny*U1;
    W=My\Uy;
    Z=(I+tauy*dt*Ay)\(-dt*W);
    U2=Z+U1;
    
    % dt/2 in x-directory
    Ux=Nx*U2;
    W=Mx\Ux;
    Z=(I+taux*dt/2*Ax)\(-dt/2*W);
    U=Z+U2;
    
end

