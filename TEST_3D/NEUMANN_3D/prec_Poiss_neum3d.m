function [U]=prec_Poiss_neum3d(f)
global K;
fc=mirt_dctn(f);
uc=fc./K;
U=mirt_idctn(uc);