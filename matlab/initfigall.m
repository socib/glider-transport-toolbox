function [] = initfigall(larg,haut,fontsiz)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fonction initialisation figure
% IN :
%   - larg, haut	: largeur et hauteur de la figure (en centim�tres)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    fontsiz = 12;
end

set(gcf,'Units','centimeters');
set(gcf,'position',[1 1 larg haut])

set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperSize',[larg haut])
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition',[0,0,1,1]);

set(gca, 'FontSize', fontsiz)

set(gcf,'color','w')

end