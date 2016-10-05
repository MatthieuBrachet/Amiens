%function show3D(timekk)
%##########################
%figure(2)
data =timekk;
data = smooth3(data,'box',5);
isoval = 0.0;
h1 = patch(isosurface(data,isoval),...
 'FaceColor','blue',...
 'EdgeColor','none',...
 'AmbientStrength',.2,...
 'SpecularStrength',.7,...
 'DiffuseStrength',.4);


isonormals(data,h1)
patch(isocaps(data,isoval),...
 'FaceColor','interp',...
 'EdgeColor','none')
colormap winter
colorbar
daspect([1,1,1])
axis tight
view(3)
camlight right
camlight left
set(gcf,'Renderer','zbuffer');
lighting phong
%saveas(gcf,f3d,'jpg')
