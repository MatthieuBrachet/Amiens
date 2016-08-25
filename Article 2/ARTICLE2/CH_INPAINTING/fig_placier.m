function [] = fig_placier

% fig_placier; (ex position.m)
% Ce programme positionne les figures presentes dans l'ecran
%
% GGB 27-03-99, 2h12
% modif 16-08-2000 : tout ecran, autant de figures qu'il y en a
% adapt. s.felix, 2004.10
% modified 30-10-04 for pbg4 with screen 1280x854

set(0,'Units','normalized')
%Screen = get(0,'ScreenSize');
xyratio=0.6672;
Screen=[0 0 xyratio 1];

% quelles figures sont ouvertes ?
if nargin<3 fenetres = sort(get(0,'children')); end

% combien de lignes, combien de colonnes pour les repartir ?
n_fig = length(fenetres);
if n_fig <= 6 nfigrow = 2;
else nfigrow = ceil(sqrt(n_fig));
end
nfigcol = ceil(n_fig/nfigrow);

% l'ordre de "fenetres" est modifie pour que les figures soit disposees "dans l'ordre"
for k = 1:nfigrow-1
fenetres((k-1)*nfigcol+1:nfigcol+(k-1)*nfigcol) = fenetres(nfigcol+(k-1)*nfigcol:-1:(k-1)*nfigcol+1);
end
fenetres((nfigrow-1)*nfigcol+1:length(fenetres)) = fenetres(length(fenetres):-1:(nfigrow-1)*nfigcol+1);

% vecteur des positions
% Corrections pour ne pas que la fenetre sorte de l'ecran faite manuellement (_o_) 
border = 0;
header = 0.09;
xoffset = 0; % offset ÌÊ partir de la gauche
yoffset = 0.0; % ofset ÌÊ partir du bas

hspace = Screen(3)/100;
vspace = Screen(4)/100;

WinWd = (Screen(3)-xoffset)/nfigcol-2*(border+hspace);
WinHt = (Screen(4)-yoffset)/nfigrow-(header+2*vspace);

for k = 1:n_fig
x = xyratio-(rem(k-1,nfigcol)+1)*(WinWd+2*border+hspace);
y = 1-((1+floor((k-1)/nfigcol))*(WinHt+2*border+vspace+header));
pos = [x y WinWd WinHt];
set(fenetres(k),'Units','normalized')
figure(fenetres(k));
set(fenetres(k),'Position',pos);
end