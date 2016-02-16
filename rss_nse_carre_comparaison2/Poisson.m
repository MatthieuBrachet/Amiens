function [Psi] = Poisson(N,W)
%test poisson - nulle au bord
% ATTENTION : on résout DDU = f et non -DDU = f!!!
h=1/(N+1);
F=-W;
[Psi] = solvpoiss2(F,N,h);
end


