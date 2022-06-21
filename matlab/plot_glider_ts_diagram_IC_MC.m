function [done] = plot_glider_ts_diagram(data,datamode,dataset,missionName,imageDir,logo,showplot) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function plot_glider_ts_diagram(data,datamode,dataset,missionName,imageDir,logo,showplot) 
%
% Purpose:                                                                           
% - Plot TS diagram for each transect
%                   
% Inputs:                                                                            
% - data        = data
% - datamode    = 'rt' o 'dt'
% - dataset     = 'glider' o 'wmop' o 'cmemsv02-med-ingv-an-fc'
% - missionName = mission name
% - imageDir    = figure output directory   
% - logo        = 'socib','copernicus','nologo'                                           
% - showplot    = 'on' or 'off' (to show plot on the screen)                                                                                                             
%   
% Outputs:
% - figures of TS diagram channels in color
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 21-Jun-2018 (Melanie Juza, mjuza@socib.es)                       
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure
if ~strcmp (logo, 'nologo')
im=get_logo([logo]);
end

% Parameters
font=17;
T0=[8:0.05:30];
S0=[34.4:0.01:38.8];
burgundy=[0.65 0.16 0.16];

% Title    
str0=['Pot. temperature / salinity diagram - ',strrep(dataset,'_','-'),' (',datamode,')'];
str1=['Balearic Channels - ', missionName];
str2=['From ' datestr(min(data.time),1) ' to '  datestr(max(data.time),1)];

% Distinguishing Ibiza and Mallorca Channels
indIC=find(data.longitude>=0.1   & data.longitude<=1.1);
indMC=find(data.longitude>=1.675 & data.longitude<=2.15);
tempIC=data.ptemp(:,indIC); saltIC=data.salinity(:,indIC);
tempMC=data.ptemp(:,indMC); saltMC=data.salinity(:,indMC);

% Min/max
minT=max(12,floor(min(min(data.ptemp))*10)/10)-0.1;
maxT=ceil(max(max(data.ptemp))*10)/10+0.1;
minS=max(36,floor(min(min(data.salinity))*10)/10);
maxS=min(38.8,ceil(max(max(data.salinity))*10)/10)+0.05;

% Figure
figure('visible',showplot); initfigall(25,20); clf;
%plot(saltIC,tempIC,'Color','k',     'Linestyle','.'); hold all 
%plot(saltMC,tempMC,'Color',burgundy,'Linestyle','.'); hold all 
plot(saltIC,tempIC,'.','Color','k'     ); hold all 
plot(saltMC,tempMC,'.','Color',burgundy); hold all
limX = [minS maxS];
limY = [minT maxT];
xlim(limX);
ylim(limY);
text(limX(1)+(limX(2)-limX(1))/15,limY(1)+(limY(2)-limY(1))/10,'Ibiza Channel','Color','k','FontSize',font);
text(limX(1)+(limX(2)-limX(1))/15,limY(1)+(limY(2)-limY(1))/20,'Mallorca Channel','Color',burgundy,'FontSize',font);
ylabel('Pot. temperature (ºC)','fontsize',font);
xlabel('Salinity','fontsize',font);
title(sprintf('%s\n%s\n%s',str0,str1,str2),'FontSize',font);
cnt_sigma(T0,S0);
set(gca,'Fontsize',font); 
grid on; 
% Logo
if ~strcmp (logo, 'nologo')
%axes('position',[0.01,0.87,0.12,0.12])
%imshow(im)
axes('position',[0.02,0.87,0.09,0.12])
image(im)  
h = gca; h.XAxis.Visible = 'off';
h = gca; h.YAxis.Visible = 'off';
end
set(gcf,'renderer','zbuffer');
% Print
imageName = ['TSdiagram_',dataset,'_',missionName,'_',datamode,'_IC_MC'];
print('-dpng',fullfile(imageDir,imageName));
display(['    - ',imageName,'.png saved in ',imageDir])

done ='ok';

end


 
 








