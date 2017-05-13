clc;
clear all;
close all;

%% *** OPTIONS ************************************************************
film='no';

%% *** STABILITY DATA *****************************************************
ddt=10^-5;
epsilon=0.01;
tau=50;
itermax=100;
lambda=75000;

%% *** TRAITEMENT DE L'IMAGE **********************************************
% lecture de l'image
image='salon.jpg';
RGB = imread(image);
N=min(size(RGB(:,:,1)));
RGB=RGB(1:N,1:N,1:3);
% calcul de l'image en n&b
for kk=1:3
    Ig = im2double(RGB(:,:,kk));
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

    %% *** ITERATIONS *********************************************************
    clear stabi iter
    stabi(1)=1; iter=0;
    while iter<itermax & stabi(end)>10^-5
        clc; kk 
        iter=iter+1
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
        end

        maxi(iter)=max(max(I));
        mini(iter)=min(min(I));
        stabi(iter)=norm(Inew-I,'inf')./norm(I,'inf');
        time(iter)=iter*ddt;

        I=Inew;
    end
    
    Im=(I+1)/2;
    Is=reshape(Im,N,N);
    RGBs(:,:,kk)=Is;
end

figure(1)
imshow(RGB)

figure(2)
imshow(RGBs)

figure(3)
for k=1:3
    subplot(2,3,k)
    imshow(RGB(:,:,k))
    
    subplot(2,3,k+3)
    imshow(RGBs(:,:,k))
end


fig_placier