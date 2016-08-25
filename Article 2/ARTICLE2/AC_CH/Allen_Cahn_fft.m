%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allen-Cahn w cosfft
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
global h N alpha beta;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %NUMERICAL DATA
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N=64;%Number of grid points in each direction 
             %N Must be odd !!!!, ideally 2^p+1
        h=1/(N-1);
        x=0:h:1;
        y=x;
        [X,Y]=meshgrid(x,y);
        
        %Time Step size
        dt=0.00001;
        nper=10;
        Nmax=500;
        epsilon=0.01;
        %Initial datum
        alpha=1;
        beta=dt;
        U=1-2*rand(N*N,1);
        
        
        k=0;
        while k < Nmax
            %computation of the r.h.s
            Fe=U.*(U.^2-1)/epsilon^2;
            F=U-dt*Fe;
            %solution of Neumann problem
            U=solv_poiss_neum(F);
%     
            if mod(k, nper)==0
                W=reshape(U,N,N);
                contourf(X,Y,W)
                drawnow
            end
            k=k+1;
        end
        
        