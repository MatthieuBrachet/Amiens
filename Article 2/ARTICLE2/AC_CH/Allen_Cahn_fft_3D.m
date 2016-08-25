%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Allen-Cahn w cosfft
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
global h N alpha beta;
global LAPLA ;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %NUMERICAL DATA
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N=32;%Number of grid points in each direction 
             %N Must be odd !!!!, ideally 2^p+1
        h=1/(N-1);
        x=0:h:1;
        y=x;
        z=x;
        [X,Y,Z]=meshgrid(x,y,z);
        


        %
        %3D Neumann Frequency Matrix
        %

  DX2=cos(pi*(0:N-1)'*ones(1,N)*h);
  DY2=cos(pi*ones(N,1)*(0:N-1)*h);
  for j=1:N
          DDX2(:,:,j)=2*(1-DX2);
          DDY2(:,:,j)=2*(1-DY2);
      end
      for k=1:N
          DDZ2(k,:,:)=2*(1-DX2);
      end
      
      LAPLA=(DDX2+DDY2+DDZ2)/h^2;
      
      
      %
      %Initial Datum
      %
      U=1-2*rand(N,N,N);
      U0=U;
      UC=dctn(U);%dct(dct(dct(U)')')';
      UC0=UC;
      
      %
      dt=0.00001;
      epsilon=0.05;
      alpha=1;beta=dt;
      K=alpha+beta*LAPLA;
      temps=[];
      Kmax=500;
      nper=10;
      k=0;
      while k< Kmax
          F=U0+dt*U0.*(1-U0.^2)/epsilon^2;
          
          %U=solv_poiss_neumann_3D(F);
          UC=dctn(F)./K;
          U=idctn(UC);
          if mod(k,nper)==0
             temps=[temps k*dt];
             figure(1)
             
             isosurface(X,Y,Z,U)
             view(3)
             figure(2)
             contourf(U(:,:,ceil((N-1)/2)))
             
           
         
             figure(5)
             UP=(U>0);
             UPP=UP.*U;
             isosurface(X,Y,Z,UPP)
             view(3)
             drawnow
             fig_placier ;
         end
           U0=U;
          k=k+1;
      end
