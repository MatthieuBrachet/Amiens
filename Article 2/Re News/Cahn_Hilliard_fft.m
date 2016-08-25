%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cahn-Hilliard w cosfft
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
global h N alpha beta;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %NUMERICAL DATA
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N=127;%Number of grid points in each direction 
             %N Must be odd !!!!, ideally 2^p+1
        h=1/(N-1);
        x=0:h:1;
        y=x;
        tau=1;
        [X,Y]=meshgrid(x,y);
       
       DF=gallery('tridiag',N);
       DF(1,1)=1;
       DF(N,N)=1;
       IDF=speye(N,N);
       LAPLA=kron(IDF,DF)+kron(DF,IDF);
       LAPLA=LAPLA/h^2;
       m=N*N;
       NDIM=m;
    
        %Time Step size
        dt=0.000001;
        nper=50;
        Nmax=5000;
        epsilon=0.01;
        % The LInear part: Block Matrix
%       MI11=speye(m,m);
%       MI12=dt*LAPLA;
%       MI21=-epsilon*LAPLA;
%       MI22=speye(m,m);
%       MM=[MI11 MI12; MI21 MI22];
      
      DX2=2*(1-cos(pi*(0:N-1)'*ones(1,N)*h))/h^2;
      DY2=2*(1-cos(pi*ones(N,1)*(0:N-1)*h))/h^2;
      KK=1+tau^2*dt*epsilon*(DX2+DY2).^2;
        
        %Initial datum
        alpha=1;
        beta=dt;
        U=1-2*rand(N*N,1);
        U=U-mean(U);
        U0=U;
        W0=U;
        W=U;
        
        ERR=[1.e-6];
        k=0;
        temps=[0];
        while k < Nmax
            %computation of the r.h.s
                 NLIN=-U.*(1-U.^2)/epsilon;

            
              alpha=1;beta=epsilon*dt;
             
              W2=epsilon*LAPLA*U0+NLIN;
              %
              W22=reshape(W2,N,N);

              
              mu2=idctn(dctn(W22)./KK);
              mu=reshape(mu2,N*N,1);
             
             %mu=bi_solv_poiss_neum(W2);
             U=U0-dt*LAPLA*mu;
             
%
            
%     
            if mod(k, nper)==0
                temps=[temps k*dt];
               
            
             
                figure(1)
                U3D=reshape(U,N,N);
                contourf(X,Y,U3D)
                
                drawnow
                fig_placier;
            end
             U0=U;
             W0=W;
            k=k+1;
        end
        
        