function [grad_I,grad_II,grad_III,grad_IV,grad_V,grad_VI]=...
    gr99(funfI,funfII,funfIII,funfIV,funfV,funfVI,n,nn)
global xi eta dxi deta
global beta betacr;
global pg kg;
global gxi_I gxi_II gxi_III gxi_IV gxi_V gxi_VI;
global geta_I geta_II geta_III geta_IV geta_V geta_VI;
global kmat
global eta_c jeta_c xi_c ixi_c
global alfasp betasp gamasp
global m1 m2
global umat1 lmat1

grad_I=zeros(nn,nn,3);grad_II=zeros(nn,nn,3);grad_III=zeros(nn,nn,3);
grad_IV=zeros(nn,nn,3);grad_V=zeros(nn,nn,3);grad_VI=zeros(nn,nn,3);
funspline=zeros(nn,1);
xc0=zeros(nn,1);xc1=zeros(nn,1);xc2=zeros(nn,1);xc3=zeros(nn,1);

%% SET 1 : 
va_fI=zeros(4*(nn-1),nn); 

% PANEL I :
va_fI(1:nn-1,1:nn)=funfI(1:nn-1,1:nn);

% PANEL II: 
for i=1:nn-1, 
    funspline(1:nn)=funfII(i,1:nn);
    sm1=(3/deta)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*deta))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*deta))*funspline(nn);
    funspline_x=zeros(nn,1);
    funspline_x(2:nn-1)=umat1\(lmat1\sm1);
    funspline_x(1)=(1/alfasp)*((betasp/deta)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*deta)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*deta))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));
    
     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/deta^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/deta;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/deta^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/deta^2;
    xc3(nn)=0;
    jj=jeta_c(i,[1:nn]);
    xeta1=eta_c(i,[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbII1=((xc3(jj).*xeta2+xc2(jj)).*xeta2+xc1(jj)).*xeta2+xc0(jj);
    va_fI(nn-1+i,1:nn)=funbII1(1:nn);
end

% PANEL III
va_fI(2*nn-1:3*nn-3,1:nn)=funfIII(1:nn-1,nn:-1:1);

%  PANEL IV:
for i=1:nn-1, 
    funspline(1:nn)=funfIV(i,1:nn);
    sm1=(3/deta)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*deta))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*deta))*funspline(nn);
    funspline_x=zeros(nn,1);
    funspline_x(2:nn-1)=umat1\(lmat1\sm1);
    funspline_x(1)=(1/alfasp)*((betasp/deta)*(funspline(2)-funspline(1))... 
        +(gamasp/(2*deta)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))...
        -(gamasp/(2*deta))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/deta^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/deta;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/deta^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/deta^2;
    xc3(nn)=0;
    jj=jeta_c(i,nn+1-[1:nn]);
    xeta1=eta_c(i,nn+1-[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbIV1=((xc3(jj).*xeta2+xc2(jj)).*xeta2+xc1(jj)).*xeta2+xc0(jj);

    va_fI(3*nn-3+i,1:nn)=funbIV1(1:nn);
end

vad_fI=(3/dxi)*pg\kg*va_fI(:,1:nn);

% for jline1=1:nn,
%     funaI=(3/dxi)*va_fI(:,jline1);
%     xsmIa=kg*funaI;
%     vad_fI(:,jline1)=pg\xsmIa;
% end

%% RESEAU 2: ASSEMBLAGE DES DONNEES SUR LE RESEAU I-BETA

vb_fI=zeros(nn,4*(nn-1));
% PANEL I :
for iline1=1:nn, 
    vb_fI(iline1,1:nn-1)=funfI(iline1,1:nn-1);
end
% FACE V:
for j=1:nn-1, 
    funspline(1:nn)=funfV(1:nn,j);
    sm1=(3/dxi)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*dxi))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*dxi))*funspline(nn);    
    funspline_x=zeros(nn,1);
    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/dxi)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*dxi)*(funspline(3)-funspline(1)))...
        - betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*dxi))*(funspline(nn-2)-funspline(nn))...
        + betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/dxi^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/dxi;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/dxi^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/dxi^2;
    xc3(nn)=0;    

    ii=ixi_c([1:nn],j);
    xxi1=xi_c([1:nn],j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaV1=((xc3(ii).*xxi2+xc2(ii)).*xxi2+xc1(ii)).*xxi2+xc0(ii);

    vb_fI(1:nn,nn-1+j)=funaV1(1:nn);
end

% FACE III:
 for iline1=1:nn,
     vb_fI(iline1,2*nn-1:3*nn-3)=funfIII(iline1,nn:-1:2); 
 end
 
% FACE VI: 
for j=1:nn-1,
    funspline(1:nn)=funfVI(1:nn,j); 
     %
    sm1=(3/dxi)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*dxi))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*dxi))*funspline(nn); 
    %
    funspline_x=zeros(nn,1);

    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/dxi)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*dxi)*(funspline(3)-funspline(1)))...
        - betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*dxi))*(funspline(nn-2)-funspline(nn))...
        + betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/dxi^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/dxi;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/dxi^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/dxi^2;
    xc3(nn)=0;
    
    ii=ixi_c(nn+1-[1:nn],j);
    xxi1=xi_c(nn+1-[1:nn],j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaVI1=((xc3(ii).*xxi2+xc2(ii)).*xxi2+xc1(ii)).*xxi2+xc0(ii);

    vb_fI(1:nn,3*nn-3+j)=funaVI1(1:nn);
end

% CALCUL DES DERIVEES BETA SUR LE RESEAU DE GRANDS CERCLES I-BETA
for iline1=1:nn,
    funbI=(3/deta)*vb_fI(iline1,:)';
    xsmIb=kg*funbI;
    vbd_fI(iline1,:)=pg\xsmIb;
end

%% SET 3: 
va_fII=zeros(4*(nn-1),nn); 
funbI1=zeros(nn,1);

% PANEL II :
for jline1=1:nn, 
    va_fII(1:nn-1,jline1)=funfII(1:nn-1,jline1);
end

% PANEL III
for i=1:nn-1, 
    funspline(1:nn)=funfIII(i,1:nn);
    %
    sm1=(3/deta)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*deta))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*deta))*funspline(nn);
    %
    funspline_x=zeros(nn,1);

    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/deta)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*deta)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*deta))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/deta^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/deta;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/deta^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/deta^2;
    xc3(nn)=0;

    jj=jeta_c(i,[1:nn]);
    xeta1=eta_c(i,[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbIII1=((xc3(jj).*xeta2+xc2(jj)).*xeta2+xc1(jj)).*xeta2+xc0(jj);

    va_fII(nn-1+i,1:nn)=funbIII1(1:nn);
end

% FACE IV: 
for jline1=1:nn,
    va_fII(2*nn-1:3*nn-3,jline1)=funfIV(1:nn-1,nn-jline1+1);    
end

% PANEL I:
for i=1:nn-1, 
   betaspline(1:nn)=beta(i,1:nn);
   funspline(1:nn)=funfI(i,1:nn);
   ppspline=spline(betaspline,funspline);
   funbI1(1:nn)=ppval(ppspline,betacr(i,nn+1-[1:nn]));
   va_fII(3*nn-3+i,1:nn)=funbI1(1:nn);
end
for i=1:nn-1,
    funspline(1:nn)=funfI(i,1:nn);
    %
    sm1=(3/deta)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*deta))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*deta))*funspline(nn);
    %
    funspline_x=zeros(nn,1);
 
    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/deta)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*deta)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*deta))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/deta^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/deta;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/deta^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/deta^2;
    xc3(nn)=0;

    jj=jeta_c(i,nn+1-[1:nn]);
    xeta1=eta_c(i,nn+1-[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbI1=((xc3(jj).*xeta2+xc2(jj)).*xeta2+xc1(jj)).*xeta2+xc0(jj);

    va_fI(3*nn-3+i,1:nn)=funbI1(1:nn);
end

for jline1=1:nn,
    funaII=(3/dxi)*va_fII(:,jline1);
    xsmIIa=kg*funaII;
    vad_fII(:,jline1)=pg\xsmIIa;
end

%% SET 4
vb_fII=zeros(nn,4*(nn-1)); 

% PANEL II : 
for iline1=1:nn,
    vb_fII(iline1,1:nn-1)=funfII(iline1,1:nn-1);
end

% PANEL V
for i=1:nn-1, 
    funspline(1:nn)=funfV(nn-i+1,1:nn);
    
    sm1=(3/dxi)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*deta))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*deta))*funspline(nn);
    
    funspline_x=zeros(nn,1);

    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/deta)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*deta)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*deta))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/deta^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/deta;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/deta^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/deta^2;
    xc3(nn)=0;    

    jj=jeta_c(i,[1:nn]);
    xeta1=eta_c(i,[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbV1=((xc3(jj).*xeta2+xc2(jj)).*xeta2+xc1(jj)).*xeta2+xc0(jj);

     vb_fII(1:nn,nn-1+i)=funbV1(1:nn);
end

% PANEL III:
for iline1=1:nn,
    vb_fII(iline1,2*nn-1:3*nn-3)=funfIV(iline1,nn:-1:2); 
end

% PANEL VI: 
for i=1:nn-1 
    funspline(1:nn)=funfVI(i,1:nn);

    sm1=(3/dxi)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*deta))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*deta))*funspline(nn);

    funspline_x=zeros(nn,1);

    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/deta)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*deta)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*deta))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/deta^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/deta;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/deta^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/deta^2;
    xc3(nn)=0;    

    jj=jeta_c(i,[1:nn]);
    xeta1=eta_c(i,[1:nn]);
    xeta2=xeta1-eta(jj);
    xeta2=xeta2';
    funbVI1=((xc3(jj).*xeta2+xc2(jj)).*xeta2+xc1(jj)).*xeta2+xc0(jj);
    %
    vb_fII(1:nn,3*nn-3+i)=funbVI1(1:nn);
end

for iline1=1:nn,
    funbII=(3/deta)*vb_fII(iline1,:)';
    xsmIIb=kg*funbII;
    vbd_fII(iline1,:)=pg\xsmIIb;
end

%% SET 5
va_fV=zeros(4*(nn-1),nn);
funaII1=zeros(nn,1);
funaIV1=zeros(nn,1);

% PANEL V
for jline1=1:nn,
    va_fV(1:nn-1,jline1)=funfV(1:nn-1,jline1);
end
% PANEL II: 
for j=1:nn-1, 
    funspline(1:nn)=funfII(1:nn,nn+1-j);
    
    sm1=(3/dxi)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*dxi))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*dxi))*funspline(nn);
    
    funspline_x=zeros(nn,1);
   
    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/dxi)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*dxi)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*dxi))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/dxi^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/dxi;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/dxi^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/dxi^2;
    xc3(nn)=0;   

    ii=ixi_c(1:nn,j);
    xxi1=xi_c(1:nn,j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaII1(1:nn)=((xc3(ii).*xxi2+xc2(ii)).*xxi2+xc1(ii)).*xxi2+xc0(ii);

    va_fV(nn-1+j,1:nn)=funaII1(1:nn);
end

% PANEL VI: 
for jline1=1:nn,
    va_fV(2*nn-1:3*nn-3,jline1)=funfVI(nn:-1:2,jline1);
end

% PANEL IV:
for j=1:nn-1,
    funspline(1:nn)=funfIV(1:nn,j);

    sm1=(3/dxi)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*dxi))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*dxi))*funspline(nn);
    
    funspline_x=zeros(nn,1);

    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/dxi)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*dxi)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*dxi))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));

     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/dxi^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/dxi;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/dxi^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/dxi^2;
    xc3(nn)=0;
    
    ii=ixi_c(1:nn,j);
    xxi1=xi_c(1:nn,j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaIV1(1:nn)=((xc3(ii).*xxi2+xc2(ii)).*xxi2+xc1(ii)).*xxi2+xc0(ii);
    
    va_fV(3*nn-3+j,1:nn)=funaIV1(1:nn);
end

for jline1=1:nn,
    funaV=(3/dxi)*va_fV(:,jline1);
    xsmVa=kg*funaV;
    vad_fV(:,jline1)=pg\xsmVa;
end

%% SET 6 :

vb_fV=zeros(nn,4*(nn-1));
funaIII1=zeros(nn,1);
funaI1=zeros(nn,1);

% PANEL V :
for iline1=1:nn, 
    vb_fV(iline1,1:nn-1)=funfV(iline1,1:nn-1);
end

% PANEL III
for j=1:nn-1, 
   funspline(1:nn)=funfIII(1:nn,nn+1-j);
   
    sm1=(3/dxi)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*dxi))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*dxi))*funspline(nn);  
   
    funspline_x=zeros(nn,1);

    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/dxi)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*dxi)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*dxi))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));
    
    % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/dxi^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/dxi;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/dxi^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/dxi^2;
    xc3(nn)=0; 
    %
    ii=ixi_c(nn+1-[1:nn],j);
    xxi1=xi_c(nn+1-[1:nn],j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaIII1(1:nn)=((xc3(ii).*xxi2+xc2(ii)).*xxi2+xc1(ii)).*xxi2+xc0(ii);
    
    vb_fV(1:nn,nn-1+j)=funaIII1(1:nn);
end

% PANEL VI :
for iline1=1:nn,
    vb_fV(iline1,2*nn-1:3*nn-3)=funfVI(nn-iline1+1,1:nn-1); 
end

% PANEL I
for j=1:nn-1, 
     funspline(1:nn)=funfI(1:nn,j);

    sm1=(3/dxi)*kmat*funspline(2:nn-1);
    sm1(1)=sm1(1)+(m1/(2*dxi))*funspline(1);
    sm1(n)=sm1(n)+(m2/(2*dxi))*funspline(nn);  

    funspline_x=zeros(nn,1);
   
    funspline_x(2:nn-1)=umat1\(lmat1\sm1);

    funspline_x(1)=(1/alfasp)*((betasp/dxi)*(funspline(2)-funspline(1))... % U_x,0
        +(gamasp/(2*dxi)*(funspline(3)-funspline(1)))...
        -betasp*funspline_x(2)-gamasp*funspline_x(3));
    funspline_x(nn)=(1/alfasp)*(-(betasp/alfasp)*(funspline(nn-1)-funspline(nn))... % U_x,N
        -(gamasp/(2*dxi))*(funspline(nn-2)-funspline(nn))...
        +betasp*funspline_x(nn-1)+gamasp*funspline_x(nn-2));   
    
     % COEFFTS C_0,C_1, C_2 AND C_3 OF CUBIC SPLINE J+1/2
    xc0(1:nn)=funspline(1:nn);
    xc1(1:nn)=funspline_x(1:nn);
    xc2(1:nn-1)=3*(funspline(2:nn)-funspline(1:nn-1))/dxi^2 ...
                -(2*funspline_x(1:nn-1)+funspline_x(2:nn))/dxi;
    xc2(nn)=0;        
    xc3(1:nn-1)=2*(funspline(1:nn-1)-funspline(2:nn))/dxi^3 ...
                + (funspline_x(1:nn-1)+funspline_x(2:nn))/dxi^2;
    xc3(nn)=0;    
    
    ii=ixi_c(nn+1-[1:nn],j);
    xxi1=xi_c(nn+1-[1:nn],j);
    xxi2=xxi1'-xi(ii);
    xxi2=xxi2';
    funaI1(1:nn)=((xc3(ii).*xxi2+xc2(ii)).*xxi2+xc1(ii)).*xxi2+xc0(ii);
    
    vb_fV(1:nn,3*nn-3+j)=funaI1(1:nn);
end

for iline1=1:nn,
    funbV=(3/deta)*vb_fV(iline1,:)';
    xsmVb=kg*funbV;
    vbd_fV(iline1,:)=pg\xsmVb;
end

%% ASSEMBLAGE
dgxi_fI=vad_fI([1:nn],[1:nn]);
dgeta_fI=vbd_fI([1:nn],[1:nn]);

dgxi_fII=vad_fII([1:nn],[1:nn]);
dgeta_fII=vbd_fII([1:nn],[1:nn]);

dgxi_fIII=vad_fI(2*(nn-1)+[1:nn],nn+1-[1:nn]);
dgeta_fIII=vbd_fI([1:nn],2*(nn-1)+nn-[1:nn]+1); 

dgxi_fIV=vad_fII(2*(nn-1)+[1:nn],nn-[1:nn]+1);
dgeta_fIV=vbd_fII([1:nn],2*(nn-1)+nn-[1:nn]+1);

dgxi_fV=vad_fV([1:nn],[1:nn]);
dgeta_fV=vbd_fV([1:nn],[1:nn]);

dgxi_fVI=vad_fV(2*(nn-1)+nn+1-[1:nn],[1:nn]);
dgeta_fVI=vbd_fV(nn-[1:nn]+1,2*(nn-1)+[1:nn]);


% 8.1-FACE I: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
grad_I(:,:,1)=dgxi_fI(:,:).*gxi_I(:,:,1) + dgeta_fI(:,:).*geta_I(:,:,1);
grad_I(:,:,2)=dgxi_fI(:,:).*gxi_I(:,:,2) + dgeta_fI(:,:).*geta_I(:,:,2);
grad_I(:,:,3)=dgxi_fI(:,:).*gxi_I(:,:,3) + dgeta_fI(:,:).*geta_I(:,:,3);

% 8.2-FACE II: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
grad_II(:,:,1)=dgxi_fII(:,:).*gxi_II(:,:,1) + dgeta_fII(:,:).*geta_II(:,:,1);
grad_II(:,:,2)=dgxi_fII(:,:).*gxi_II(:,:,2) + dgeta_fII(:,:).*geta_II(:,:,2);
grad_II(:,:,3)=dgxi_fII(:,:).*gxi_II(:,:,3) + dgeta_fII(:,:).*geta_II(:,:,3);

% 8.3-FACE III: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA 
grad_III(:,:,1)=dgxi_fIII(:,:).*gxi_III(:,:,1) - dgeta_fIII(:,:).*geta_III(:,:,1);
grad_III(:,:,2)=dgxi_fIII(:,:).*gxi_III(:,:,2) - dgeta_fIII(:,:).*geta_III(:,:,2);
grad_III(:,:,3)=dgxi_fIII(:,:).*gxi_III(:,:,3) - dgeta_fIII(:,:).*geta_III(:,:,3);

% 8.4 - FACE IV: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
grad_IV(:,:,1)=dgxi_fIV(:,:).*gxi_IV(:,:,1) - dgeta_fIV(:,:).*geta_IV(:,:,1);
grad_IV(:,:,2)=dgxi_fIV(:,:).*gxi_IV(:,:,2) - dgeta_fIV(:,:).*geta_IV(:,:,2);
grad_IV(:,:,3)=dgxi_fIV(:,:).*gxi_IV(:,:,3) - dgeta_fIV(:,:).*geta_IV(:,:,3);

% 8.5 - FACE V: CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
grad_V(:,:,1)=dgxi_fV(:,:).*gxi_V(:,:,1) + dgeta_fV(:,:).*geta_V(:,:,1);
grad_V(:,:,2)=dgxi_fV(:,:).*gxi_V(:,:,2) + dgeta_fV(:,:).*geta_V(:,:,2);
grad_V(:,:,3)=dgxi_fV(:,:).*gxi_V(:,:,3) + dgeta_fV(:,:).*geta_V(:,:,3);

% 8.6 - FACE VI : CALCUL DU GRADIENT EN FONCTION DES DERIVEES ALPHA ET BETA
grad_VI(:,:,1)=-dgxi_fVI(:,:).*gxi_VI(:,:,1) + dgeta_fVI(:,:).*geta_VI(:,:,1);
grad_VI(:,:,2)=-dgxi_fVI(:,:).*gxi_VI(:,:,2) + dgeta_fVI(:,:).*geta_VI(:,:,2);
grad_VI(:,:,3)=-dgxi_fVI(:,:).*gxi_VI(:,:,3) + dgeta_fVI(:,:).*geta_VI(:,:,3);

%% 1/2 SOMME ARRETES, 1/3 SOMME SOMMETS
% COMPONENT 1
uwk_I=grad_I(1:nn,1:nn,1);uwk_II=grad_II(1:nn,1:nn,1);uwk_III=grad_III(1:nn,1:nn,1);
uwk_IV=grad_IV(1:nn,1:nn,1);uwk_V=grad_V(1:nn,1:nn,1);uwk_VI=grad_VI(1:nn,1:nn,1);
[grad_I(1:nn,1:nn,1),grad_II(1:nn,1:nn,1),grad_III(1:nn,1:nn,1),grad_IV(1:nn,1:nn,1),grad_V(1:nn,1:nn,1),grad_VI(1:nn,1:nn,1)]=...
    ds98(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
% COMPONENT 2
uwk_I=grad_I(1:nn,1:nn,2);uwk_II=grad_II(1:nn,1:nn,2);uwk_III=grad_III(1:nn,1:nn,2);
uwk_IV=grad_IV(1:nn,1:nn,2);uwk_V=grad_V(1:nn,1:nn,2);uwk_VI=grad_VI(1:nn,1:nn,2);
[grad_I(1:nn,1:nn,2),grad_II(1:nn,1:nn,2),grad_III(1:nn,1:nn,2),grad_IV(1:nn,1:nn,2),grad_V(1:nn,1:nn,2),grad_VI(1:nn,1:nn,2)]=...
    ds98(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);
% COMPONENT 3
uwk_I=grad_I(1:nn,1:nn,3);uwk_II=grad_II(1:nn,1:nn,3);uwk_III=grad_III(1:nn,1:nn,3);
uwk_IV=grad_IV(1:nn,1:nn,3);uwk_V=grad_V(1:nn,1:nn,3);uwk_VI=grad_VI(1:nn,1:nn,3);
[grad_I(1:nn,1:nn,3),grad_II(1:nn,1:nn,3),grad_III(1:nn,1:nn,3),grad_IV(1:nn,1:nn,3),grad_V(1:nn,1:nn,3),grad_VI(1:nn,1:nn,3)]=...
    ds98(uwk_I,uwk_II,uwk_III,uwk_IV,uwk_V,uwk_VI,n,nn);