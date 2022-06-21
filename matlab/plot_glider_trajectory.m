function [done] = plot_glider_trajectory(lon,lat,missionName,datamode,imageDir,logo,showplot)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function plot_glider_trajectory(lon,lat,missionName,datamode,imageDir,logo,showplot)  
%                                                                     
% Purpose:                                                                           
% - Plot glider trajectory during the mission
%                   
% Inputs:                                                                            
% - longitude 
% - latitude
% - missionName = mission name
% - datamode    = 'rt' o 'dt'                                                           
% - imageDir    = figure output directory     
% - logo       = 'socib','copernicus','nologo'
% - showplot  = 'on' or 'off' (to show plot on the screen)                                                                                                             
%   
% Outputs:
% - figure of glider trajectory
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 19-June-2018 (Melanie Juza, mjuza@socib.es)                       
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure
font=17;
if ~strcmp (logo, 'nologo')
im=get_logo([logo]);
end
strtitle=['All transects ',missionName,' (',datamode,')'];

%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('visible',showplot); initfigall(25,20); clf;
m_proj('mercator','longitudes',[0 3],'latitudes',[38 40]);
m_usercoast('medcoast_huge_noriv','patch',[.4 .4 .4]);
m_usercoast('medcoast_huge_noriv','linewidth',1,'color','k');
m_grid2('box','on','tickdir','in','fontsize', font);
hold on;
m_plot( lon, lat,'.b');
title (strtitle,'fontsize', font)
set(gca,'fontsize', font)
% Logo
if ~strcmp (logo, 'nologo')
 if strcmp (logo, 'copernicus')
  axes('position',[0.01,0.93,0.085,0.085])
 else
  %axes('position',[0.005,0.92,0.07,0.07])
  axes('position',[0.01,0.92,0.055,0.07])
 end
%imshow(im)
image(im)  
h = gca; h.XAxis.Visible = 'off';
h = gca; h.YAxis.Visible = 'off';
end
set(gcf,'renderer','zbuffer');
% Print
imageName = ['mapAllTransects_',missionName,'_',datamode];
print('-dpng',fullfile(imageDir,imageName)); 
display(['    - ',imageName,'.png saved in ',imageDir])
done = i;

end

