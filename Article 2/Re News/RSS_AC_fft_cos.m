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
        N=128;%Number of grid points in each direction 
             %N Must be odd !!!!, ideally 2^p+1
        h=1/(N-1);
        x=0:h:1;
        y=x;
        [X,Y]=meshgrid(x,y);
        %
        %Building of teh compact schem matrix
        %
        [P,Q]=CS_lapla_matrix_Neumann_bc(N-2,h);
        LP4=P\Q;
        ID=speye(N,N);
        LAPLA4=kron(LP4,ID)+kron(ID,LP4);
        %Time Step size
        dt=0.00001;
        nper=10;
        Nmax=2000;
        epsilon=0.01;
        %Initial datum
         tau=1;
        alpha=1;
        beta= tau*dt;
       
        U=1-2*rand(N*N,1);
        
        
        k=0;
        while k < Nmax
            %computation of the r.h.s
            Fe=U.*(U.^2-1)/epsilon^2;
            F=-dt*LAPLA4*U-dt*Fe;
            %solution of Neumann problem
            DELTA=solv_poiss_neum(F);
            U=U+DELTA;
%     
            if mod(k, nper)==0
                W=reshape(U,N,N);
                contourf(X,Y,W)
                drawnow
            end
            k=k+1;
        end
        
        