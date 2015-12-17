function [bwx,bwy] = bord(psi)
%mise à jour des conditions de bord
global h Kx Ky bb
bwx=Kx*psi;
bwy=Ky*psi+25/(6*h)*bb;
end

