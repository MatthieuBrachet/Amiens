function  [u_val] = solvpoiss(F,n,dt,h,tau,Re)

% solveur de :
%      (Id + dt*tau/Re A2) x = F
%

nx=n;
ny=n;
x=1:nx;
y=1:ny;
dx =h; dy=h;
F=reshape(F,nx,ny);
%on passe dans la base des vecteurs propres du laplacien
%via une FFT 2D de sinus
u = dst2(F);    
%valeurs propres de la matrice du laplacien
d = 4*(sin(x'/2*pi/(nx+1)).^2*ones(1,ny)/dx^2+ ...
                 ones(nx,1)*sin(y/2*pi/(ny+1)).^2/dy^2);
d=1+(tau*dt/Re)*d;
%resolution dna sla base des vecteurs propres              
u = u./d;   
%reecriture du resultat dans la base canononque via FFT inverse
u = dst2(u); 
        
 
u_val=u(:);
return