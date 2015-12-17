%% sous figure

clc
close all

nb_contour=1000;

[ PP1, PP2, PP3, PP4, PP5, PP6, PP7, PP8, PP9 ] = vortex(PP);
[ X1, X2, X3, X4, X5, X6, X7, X8, X9 ] = vortex(X);
[ Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9 ] = vortex(Y);


figure(1)

subplot(331)
contour(X7,Y7,PP7,nb_contour)
title('zone 7')

disp('---------ZONE 7---------')
[MIN,MAX] = MinMax(PP7,X7,Y7);
disp('max : ')
MAX
disp('min : ')
MIN


subplot(332)
contour(X8,Y8,PP8,nb_contour)
title('zone 8')

disp('---------ZONE 8---------')
[MIN,MAX] = MinMax(PP8,X8,Y8);
disp('max : ')
MAX
disp('min : ')
MIN

subplot(333)
contour(X9,Y9,PP9,nb_contour)
title('zone 9')

disp('---------ZONE 9---------')
[MIN,MAX] = MinMax(PP9,X9,Y9);
disp('max : ')
MAX
disp('min : ')
MIN

subplot(334)
contour(X4,Y4,PP4,nb_contour)
title('zone 4')

disp('---------ZONE 4---------')
[MIN,MAX] = MinMax(PP4,X4,Y4);
disp('max : ')
MAX
disp('min : ')
MIN

subplot(335)
contour(X5,Y5,PP5,nb_contour)
title('zone 5')

disp('---------ZONE 5---------')
[MIN,MAX] = MinMax(PP5,X5,Y5);
disp('max : ')
MAX
disp('min : ')
MIN

subplot(336)
contour(X6,Y6,PP6,nb_contour)
title('zone 6')

disp('---------ZONE 6---------')
[MIN,MAX] = MinMax(PP6,X6,Y6);
disp('max : ')
MAX
disp('min : ')
MIN


subplot(337)
contour(X1,Y1,PP1,nb_contour)
title('zone 1')

disp('---------ZONE 1---------')
[MIN,MAX] = MinMax(PP1,X1,Y1);
disp('max : ')
MAX
disp('min : ')
MIN

subplot(338)
contour(X2,Y2,PP2,nb_contour)
title('zone 2')

disp('---------ZONE 2---------')
[MIN,MAX] = MinMax(PP2,X2,Y2);
disp('max : ')
MAX
disp('min : ')
MIN

subplot(339)
contour(X3,Y3,PP3,nb_contour)
title('zone 3')

disp('---------ZONE 3---------')
[MIN,MAX] = MinMax(PP3,X3,Y3);
disp('max : ')
MAX
disp('min : ')
MIN


