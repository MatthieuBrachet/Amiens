%
%test
%
clear all;close all;
global N h;
%
N=64;
h=1/(N-1);
alpha=1;
beta=1;
%
x=0:h:1;
[X,Y,Z]=meshgrid(x,x,x);
%
%% 3D Neumann Frequency Matrix
% one=ones(N,1);
% K=cos(pi*(0:N-1)*h);
% %
% DX2=kron(one,K);
% DY2=DX2';
% 
% 
% for j=1:N
%      DDX2(:,:,j)=2*(1-DX2);
%      DDY2(:,:,j)=2*(1-DY2);
%      DDZ2(j,:,:)=2*(1-DX2);
% end
%       
% LAPLA=(DDX2+DDY2+DDZ2)/h^2;
% K=alpha+beta*LAPLA.^2;

%%%%%%%%%%%%%%%%%%%%%
uex=cos(pi*X).*cos(pi*Y).*cos(pi*Z);
Uex=reshape(uex,N*N*N,1);

f=(alpha+beta*9*pi^4).*uex;
F=reshape(f,N*N*N,1);
%%%%%%%%%%%%%%%%%%%%%
%R�solustion du pb de Poisson, sch�mas comapcts
%% solver classique

% [a0,P,Q] = Mlaplacien2(N-2,2);
% A=P\Q;
% id=speye(N,N);
% Ax=kron(kron(A,id),id);
% Ay=kron(kron(id,A),id);
% Az=kron(kron(id,id),A);
% At=sparse(Ax+Ay+Az);
% ID=speye(size(At));
% MAT=alpha*ID+beta*At^2;
% 
% Uclassic=MAT\F;

%% solver precond.
uprec=prec_bilap_neum3d(F,alpha,beta);
Uprec=reshape(uprec,N*N*N,1);

%% results
disp(' *** Error ***** ')
disp([' when N = ',num2str(N)])
disp([' precond. :' num2str(norm(Uprec-Uex,'inf'))])
%disp([' classic. :' num2str(norm(Uclassic-Uex,'inf'))])
disp(' ************** ')
