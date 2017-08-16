clc;
clear all;
close all;

%% *** OPTIONS ************************************************************
film='no';
image='debussy1'

%% *** STABILITY DATA *****************************************************
ddt=3*10^-4;
epsilon=0.03;
tau=1.0;
itermax=1000;
nper=1;
%lambda=750000;
lambda=10^10;
%% *** TRAITEMENT DE L'IMAGE **********************************************
% lecture de l'image
RGB = imread([image '.jpeg']);
N=min(size(RGB(:,:,1)));
h=1/(N-1);
RGB=RGB(1:N,1:N,1:3);
% calcul de l'image en n&b
Ig = im2double(rgb2gray(RGB));
m=min(min(Ig)); M=max(max(Ig));
f0=(Ig-m)./(M-m);
figure(1000)
imshow((2*rand(N,N)-1).*f0)
title(['segmented image ','t=',num2str(0,'%.8f')]);
pause
f0sq=f0;
f0=reshape(f0,[],1);
%I=2.*f0-1;
I=reshape((2*rand(N,N)-1).*f0sq,[],1);

%% *** SPACE DATA *********************************************************
[~,id,AA] = Mlaplacien2(N-2,2);
Ay=kron(AA,id);
Ax=kron(id,AA);
A=Ax+Ay;

[~,MM,NN] = Mlaplacien2(N-2,4);
My=sparse(kron(MM,id));
Mx=sparse(kron(id,MM));
%P=Mx+My;

Ny=sparse(kron(NN,id));
Nx=sparse(kron(id,NN));
%Q=Nx+Ny;

%A4=Mx\Nx+My\Ny;%P\Q;
%INTERGRATOR TRAPEZOIDAL RULE MATRIX
TR=ones(N,N);
TR(1,1)=1/4;TR(1,2:N-1)=1/2;TR(1,N)=1/4;
TR(2:N-1,1)=1/2;TR(2:N-1,N)=1/2;
TR(N,1)=1/4;TR(N,2:N-1)=1/2;TR(N,N)=1/4;

ID=speye(size(Ax));
%
%FFT ARRAYS
%
      %DX2=2*(1-cos(pi*(0:N-1)'*ones(1,N)*h))/h^2;
      %DY2=2*(1-cos(pi*ones(N,1)*(0:N-1)*h))/h^2;
      %KK=1+ddt*tau*(DX2+DY2);
      p = 0:N-1;
      q = 0:N-1;
      [p,q]=meshgrid(p,q);
      KK=1+ddt*(4-(2 * cos(pi * p /(N-1)) + 2 * cos(pi * q / (N-1))))/h^2;
      %
   if strcmp(film,'yes')==1
    vidObj=VideoWriter([image '.avi']);
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end


%% *** ITERATIONS *********************************************************
stabi(1)=1; iter=0;
ccc=exp(-2*ddt/(epsilon^2));
while iter<itermax & stabi(end)>10^-4
%      if stabi(end)<10^-3
%          %ddt=1.e-4;
%          epsilon=10;
%         %ddt=1.1*ddt;
%          ccc=exp(-2*ddt/(epsilon^2));
%      end
    clc; disp([iter stabi(end)]);
    iter=iter+1;
     % *** Segmentation part
     Image=reshape(I,N,N);
     Ii0=sum(TR.*(f0sq.*(1+Image)));
     Ii1=sum(TR.*(1+Image));
     Ii2=sum(TR.*(f0sq.*(1-Image)));
     Ii3=sum(TR.*(1-Image));
     c1=Ii0/Ii1;
     c2=Ii2/Ii3;
    %c1=sum(f0.*(1+I))./sum(1+I);
    %c2=sum(f0.*(1-I))./sum(1-I);
    %
    GV=sparse(((f0-c1).^2+(f0-c2).^2));
    G=ID+ddt*lambda*diag(GV);
    I1=G\(I-lambda*ddt*((f0-c1).^2-(f0-c2).^2));
    %I1=I-ddt*lambda*((1+I).*(f0-c1).^2-(1-I).*(f0-c2).^2);
    % *** Heat equation part
   % Z=(ID+tau*ddt*A)\(I1-ddt*A*I1);
   %Computation of A4*I1
   I10=Mx\(Nx*I1);
   I11=My\(Ny*I1);
   ZZ=I10+I11;
   ZZ=-ddt*A*I1;
   %
   %Resolution par fft_cos
   %
    Z2=reshape(ZZ,N,N); 
    %Z3=dct2(Z2);
    %Z3(1,1)=0;
              mu2=idct2(dct2(Z2)./KK);
              Z=reshape(mu2,N*N,1);
              
    %Z=(ID+tau*ddt*A)\(-ddt*ZZ);
    I2=I1+Z;
  %  I2=(ID+ddt*A)\(I1);
    %I2=Z+I1;
    % *** Allen-Cahn non linear part
    
    num=I2;
    
    denom=sqrt(ccc+I2.^2.*(1-ccc));
    Inew=num./denom;

    

   
    
   % G=-lambda*((1+I3).*(f0-c1).^2-(1-I3).*(f0-c2).^2);
   % Inew=I3-ddt*G;
    
   if mod(iter,nper)==0
       niv1=max(max(Inew));
       niv2=min(min(Inew));
       niv=(niv1+niv2)/2;
       Imi=(Inew<niv);
Ipl=(Inew>niv);
%post-traitement
%Im=(I-(Imi+Ipl)/2)./(Ipl-Imi);
Im=Ipl-Imi;

Is=reshape(Inew,N,N);

figure(1000)
imshow(Is)
title(['segmented image ','t=',num2str(iter*ddt,'%.8f')]);

drawnow
   end
   
    if strcmp(film,'yes') == 1
        Im=(Inew+1)/2;
        Is=reshape(Im,N,N);
        imshow(Is)
        title(['iteration : ' num2str(iter)]); 
        
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
    end
    
    maxi(iter)=max(max(I));
    mini(iter)=min(min(I));
    stabi(iter)=norm(Inew-I,'inf')./norm(I,'inf');
    time(iter)=iter*ddt;
    
    I=Inew;
end

if strcmp(film,'yes') == 1
    close(vidObj);
end

%% *** IMAGES *************************************************************

Imi=(I<0);
Ipl=(I>0);
%post-traitement
%Im=(I-(Imi+Ipl)/2)./(Ipl-Imi);
Im=Ipl-Imi;

Is=reshape(Im,N,N);

figure(1)
imshow(Is)
title('segmented image')
   
figure(2)
imshow(RGB)
title('original image')

figure(3)
imshow(Ig)
title('image b&w')

figure(4)
plot(time/ddt, maxi, time/ddt, mini)
legend('max','min')

figure(5)
semilogy(time/ddt, stabi)
figure(6)
imshow(-Is)
title('inverse segmented image')
fig_placier