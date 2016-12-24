% video
[n1,n2,n3,n4]=size(im_v);

% Prepare the new file.
vidObj = VideoWriter('peaks.avi');
open(vidObj);

% Create an animation.
im=squeeze(im_v(:,:,:,1));
isosurface(X,Y,Z,im);
axis([0 1 0 1 0 1])
grid on

set(gca,'nextplot','replacechildren');


for iter = 1:n4
    clc; clf; disp([num2str(iter) ' on ' num2str(n4)])
    im=squeeze(im_v(:,:,:,iter));
    im=triche(im);
    isosurface(X,Y,Z,im);
    axis([0 1 0 1 0 1])
    grid on
    hold off

    % Write each frame to the file.
    currFrame = getframe;
    writeVideo(vidObj,currFrame);
end

% Close the file.
close(vidObj);