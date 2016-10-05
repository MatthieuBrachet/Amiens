clc; clear all; close all;
% CNRSS !

% si film = 1 : faire le film,
%    film = 0 : ne pas faire de film.
film = 0;

ref=floor(10000*now);
if film==1
    % options de film
    %nbim=1000;%(=5 when dt=0.01)
    
    mkdir(['./video-' date ])
    mov=avifile(['./video-' date '/ref_' num2str(ref) '_AC.avi'],'compression','None');
    fig=figure;
end

Tmax=0.4;
ddt=10^-3;

N=20;
h=1/(N+1);
x=[0:h:1]';
[X,Y,Z]=meshgrid(x,x,x);

%U=2*rand(size(X))-1;
U=cos(pi*X).*cos(pi*Y).*cos(pi*Z);
U=reshape(U,[],1);
epsilon=0.01;

%% schéma d'ordre 4
[a0,P,Q] = Mlaplacien2(N,4);
id=speye(N+2,N+2);

Px=kron(kron(P,id), id);
Py=kron(kron(id,P), id);
Pz=kron(kron(id,id), P);

Qx=kron(kron(Q,id), id);
Qy=kron(kron(id,Q), id);
Qz=kron(kron(id,id), Q);

Zx=Px\Qx;
Zy=Py\Qy;
Zz=Pz\Qz;
%% schéma d'ordre 2
[a0,P,Q] = Mlaplacien2(N,2);
A=P\Q;

Ax=kron(kron(A,id), id);
Ay=kron(kron(id,A), id);
Az=kron(kron(id,id), A);
Id=speye(size(Ax));

t=0;
while t<Tmax
    clc; t=t+ddt
    
    % step 1
    W=(Id+ddt*Ax)\(-ddt*(Zx*U));
    V1=W+U;
    
    % step 2
    W=(Id+ddt*Ay)\(-ddt*(Zy*V1));
    V2=W+V1;
    
    % step 3
    W=(Id+ddt*Az)\(-ddt*Zz*V2);
    V3=W+V2;
    
    % step 4
    nom=V3;
    denom=V3.^2+(1-V3.^2).*exp(-2*ddt/(epsilon^2));
    denom=sqrt(denom);
    U=nom./denom;
    
    if film == 1
        % plot
        figure(1)
        isosurface(X,Y,Z,reshape(U,size(X)))
        title(['time : ', num2str(t)], 'Units', 'normalized','Position', [1 1], 'HorizontalAlignment', 'right')
        %title(['solution at time : ', num2str(t)], 'HorizontalAlignment', 'right');
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis([0 1 0 1 0 1])
    
        frame = getframe(fig);
        mov = addframe(mov,frame); 
        hold off
        
        close (1)
    end

end

if film==1
    %close(fig)
    mov=close(mov);
    
    data = fopen('AAA_VIDEO_SAVE.txt','a');
    fprintf(data,'%s\n',['date : ', date]);
    fprintf(data,'%s\n',['ref. : ', num2str(ref)]);
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','---------- numerical data ---------');
    fprintf(data,'%s\n',['grid              : ', num2str(N)] );
    fprintf(data,'%s\n',['time step         : ', num2str(ddt)] );
    fprintf(data,'%s\n','-------- mathematical data --------');
    fprintf(data,'%s\n',['epsilon           : ', num2str(epsilon)] );
    fprintf(data,'%s\n','------------ results --------------');
    fprintf(data,'%s\n',['final time        : ', num2str(t)] );
    fprintf(data,'%s\n','***********************************');
    fprintf(data,'%s\n','  ');
    fprintf(data,'%s\n','  ');
    fclose(data);
end

figure(1)
isosurface(X,Y,Z,reshape(U,size(X)))