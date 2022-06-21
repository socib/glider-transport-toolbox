function [done] = plot_glider_density_individual_sections(densMatrix,gvTimeTrans,depthLevels,lonGrid,gvbathy,missionName,transectNo,imageDir,section,datamode,dataset,logo,showplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function plot_glider_density_individual_sections(densMatrix,gvTimeTrans,depthLevels,lonGrid,gvbathy,missionName,transectNo,imageDir,section,datamode,dataset,logo,showplot)
%                                                                     
% Purpose:                                                                           
% - Plot vertical section density for each transect
%                   
% Inputs:                                                                            
% - densMatrix  = density 
% - gvTimeTrans = time  
% - depthLevels = depth   
% - lonGrid     = longitudes
% - bathy       = bathymetry
% - missionName = mission name
% - transectNo  = transect number identification
% - imageDir    = figure output directory   
% - section     = IbizaChannel  
% - datamode    = 'rt' o 'dt'
% - dataset     = 'glider' o 'wmop' o 'cmemsv02-med-ingv-an-fc'
% - logo      = 'socib','copernicus','nologo'                                           
% - showplot    = 'on' or 'off' (to show plot on the screen)                                                                                                             
%   
% Outputs:
% - figures of vertical sections of density with bathymetry line 
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 23-Jul-2018 (Melanie Juza, mjuza@socib.es)                       
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp (logo, 'nologo')
im=get_logo([logo]);
end

if strcmp (section, 'IbizaChannel')
 lonmin=0.1;
 lonmax=1.1;
elseif strcmp (section,'MallorcaChannel')
 lonmin=1.675;
 lonmax=2.15;
end

font=17;

%%% Figures of density

% Title
str0=['Density - ',strrep(dataset,'_','-'),' (',datamode,')'];
str1=[section,' - ', missionName, ' - T',num2str(transectNo)];
str2=['From ' datestr(min(gvTimeTrans),1) ' to '  datestr(max(gvTimeTrans),1)];

% Limits
minV=min(densMatrix(:));
maxV=max(densMatrix(:));

% Figure
figure('visible',showplot); initfigall(30,20); clf;
colormap('jet');
%contourGrid = ones(size(depthLevels))'*lonGrid;
%pressGrid = depthLevels(:) * ones(size(lonGrid));
%sigma = densMatrix(:,:) - 1000;
%[C,h] = contour(contourGrid,pressGrid,sigma(:,:),...
%[28.00 28.40 28.80  28.90  29.00 29.05 ],'-k', 'linewidth', 1);  hold on
%clabel(C,h,'Fontsize',12,'Color','k');hold on
pcolor(lonGrid,depthLevels,densMatrix);shading interp;
hold on;plot(lonGrid,gvbathy,'Color','k','linewidth',1);
fill( [lonGrid fliplr(lonGrid)],  [gvbathy fliplr(1100*ones(1,length(gvbathy)))], [.8 .8 .8]);
axis ij
hc= colorbar;
set(get(hc, 'xLabel'), 'String','kg/m3','Fontsize',font);    
xlabel('Longitude','Fontsize',font);
ylabel('Depth (m)','Fontsize',font);
axis([lonmin lonmax 0 1000 ])
title(sprintf('%s\n%s\n%s',str0,str1,str2),'FontSize',font);
caxis([minV maxV]); 
set(gca,'Fontsize',font); 
grid on; 
% Contour sigma
contourGrid = ones(size(depthLevels))'*lonGrid;
pressGrid = depthLevels(:) * ones(size(lonGrid));
sigma = densMatrix(:,:) - 1000;
[C,h] = contour(contourGrid,pressGrid,sigma(:,:),...
[28.00 28.40 28.80  28.90  29.00 29.05 ],'-k', 'linewidth', 1);  hold on
clabel(C,h,'Fontsize',12,'Color','k');hold on
% Logo
if ~strcmp (logo, 'nologo')
  if strcmp (logo, 'copernicus') 
    axes('position',[0.01,0.9,0.1,0.1])
  elseif strcmp (logo, 'socib') 
    %axes('position',[0.005,0.87,0.12,0.12])
    axes('position',[0.015,0.87,0.075,0.12])
  end
  %imshow(im)
  image(im)  
  h = gca; h.XAxis.Visible = 'off';
  h = gca; h.YAxis.Visible = 'off';
end
set(gcf,'renderer','zbuffer');
% Print
imageName=['Density_section_',dataset,'_',section,'_',missionName,'_T',num2str(transectNo),'_',datamode];
print('-dpng',fullfile(imageDir,imageName));
display(['    - ',imageName,'.png saved in ',imageDir])

done ='ok';

end
  
