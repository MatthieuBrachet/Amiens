clc;
clear all;
close all;

%% *** OPTIONS ************************************************************
film='no';
image='sat2';

%% *** STABILITY DATA *****************************************************
ddt=10^-5;
epsilon=0.03;
tau=50;
itermax=50;
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
Ax=kron(AA,id);
Ay=kron(id,AA);
A=Ax+Ay;

[~,MM,NN] = Mlaplacien2(N-2,2);
Mx=kron(MM,id);
My=kron(id,MM);
P=Mx+My;

Nx=kron(NN,id);
Ny=kron(id,NN);
Q=Nx+Ny;

A4=P\Q;

ID=speye(size(Ax));

if strcmp(film,'yes')==1
    vidObj=VideoWriter([image '.avi']);
    open(vidObj);
    set(gca,'nextplot','replacechildren');
end

%% *** ITERATIONS *********************************************************
stabi(1)=1; iter=0;
while iter<itermax & stabi(end)>10^-4
    clc; disp([iter stabi(end)]);
    iter=iter+1;
    % *** Allen-Cahn non linear part
    num=I;
    denom=sqrt(exp(-2*ddt/(epsilon^2))+I.^2.*(1-exp(-2*ddt/(epsilon^2))));
    I1=num./denom;

    % *** Heat equation part
    Z=(ID+tau*ddt*A)\(-ddt*A4*I1);
    I3=I1+Z;

    % *** Segmentation part
    c1=sum(f0.*(1+I3))./sum(1+I3);
    c2=sum(f0.*(1-I3))./sum(1-I3);
    
    G=-lambda*((1+I3).*(f0-c1).^2-(1-I3).*(f0-c2).^2);
    Inew=I3-ddt*G;
    
    if strcmp(film,'yes') == 1
        Im=(I+1)/2;
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
Im=(I+1)/2;
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