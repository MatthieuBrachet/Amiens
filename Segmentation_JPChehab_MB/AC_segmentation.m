clc;
clear all;
close all;

%% *** OPTIONS ************************************************************
film='no';
image='europa_one';

%% *** STABILITY DATA *****************************************************
ddt=10^-5;
epsilon=0.04;
tau=1.;
itermax=1000;
nper=1;
lambda=75000;

%% *** TRAITEMENT DE L'IMAGE **********************************************
% lecture de l'image
RGB = imread([image '.jpg']);
N=min(size(RGB(:,:,1)));
RGB=RGB(1:N,1:N,1:3);
% calcul de l'image en n&b
Ig = im2double(rgb2gray(RGB));
m=min(min(Ig)); M=max(max(Ig));
f0=(Ig-m)./(M-m);
f0=reshape(f0,[],1);
I=2.*f0-1;


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

ID=speye(size(Ax));

if strcmp(film,'yes')==1
    vidObj=VideoWriter([image '.avi']);
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end

%% *** ITERATIONS *********************************************************
stabi(1)=1; iter=0;
ccc=exp(-2*ddt/(epsilon^2));
while iter<itermax & stabi(end)>10^-4
    clc; disp([iter stabi(end)]);
    iter=iter+1;
     % *** Segmentation part
    c1=sum(f0.*(1+I))./sum(1+I);
    c2=sum(f0.*(1-I))./sum(1-I);
    %
    GV=sparse(((f0-c1).^2+(f0-c2).^2));
    G=ID+ddt*lambda*diag(GV);
    I1=G\(I-lambda*ddt*((f0-c1).^2-(f0-c2).^2));
    %I1=I-ddt*lambda*((1+I).*(f0-c1).^2-(1-I).*(f0-c2).^2);
    % *** Heat equation part
   % Z=(ID+tau*ddt*A)\(I1-ddt*A*I1);
   %Computation of A4*I1
   %I10=Mx\(Nx*I1);
   %I11=My\(Ny*I1);
   %ZZ=I10+I11;
   ZZ=A*I1;
    Z=(ID+tau*ddt*A)\(-ddt*ZZ);
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
title('segmented image')
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

fig_placier