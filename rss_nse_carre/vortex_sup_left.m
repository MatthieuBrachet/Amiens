% vortex supérieur gauche et inferieur droit

clc
close all
format shorte

[ PP1, PP2, PP3, PP4, PP5, PP6, PP7, PP8, PP9 ] = vortex(PP);
[ X1, X2, X3, X4, X5, X6, X7, X8, X9 ] = vortex(X);
[ Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9 ] = vortex(Y);

%**************************************************************************

hight_left=level(PP7,-.004,1);

figure(1)
hold on
[CC,hh]=contour(X7,Y7,PP7,20);
clabel(CC,hh);
[CC,hh]=contour(X7,Y7,hight_left,10);
clabel(CC,hh);
hold off
xlabel('x')
ylabel('y')
title('vortex superieur gauche')


[MIN,MAX] = MinMax(hight_left,X7,Y7);

disp('---- vortex superieur gauche ----')
MIN
MAX

%**************************************************************************

bottom_right=level(PP3,2.5028*10^-8,1);

figure(2)
hold on
[CC,hh]=contour(X3,Y3,PP3,25);
clabel(CC,hh);
[CC,hh]=contour(X3,Y3,bottom_right,10);
clabel(CC,hh);
hold off
xlabel('x')
ylabel('y')
title('vortex inferieur droit')


[MIN,MAX] = MinMax(bottom_right,X3,Y3);

disp('---- vortex inferieur droit ----')
MIN
MAX

%**************************************************************************

center=PP;

figure(3)
hold on
[CC,hh]=contour(X,Y,center,25);
clabel(CC,hh);
hold off
xlabel('x')
ylabel('y')
title('vortex central')


[MIN,MAX] = MinMax(center,X,Y);

disp('---- vortex central ----')
MIN
MAX





