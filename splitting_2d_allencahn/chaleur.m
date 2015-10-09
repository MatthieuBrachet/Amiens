function [ U ] = chaleur( U , dt)
    global A2 I P Q tau
    
    Wtemp=P\U;
    W=Q*Wtemp;
    Z=(I+tau*dt*A2)\(-dt*W);
    U=Z+U;

end

