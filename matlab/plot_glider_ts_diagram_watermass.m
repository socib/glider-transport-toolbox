function plot_glider_ts_diagram_watermass(gvPtemp,gvSalt,gvTimeTrans,AWo,AWr,LIW,WIW,WDW,critWIW,missionName,transectNo, imageDir, section,datamode,dataset,logo,showplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function plot_glider_ts_diagram_watermass(gvPtemp,gvSalt,AWo,AWr,LIW,WIW,WDW,critWIW,missionName,transectNo, imageDir, section,datamode,dataset,logo,showplot)
%
% Purpose:                                                                           
% - Plot TS diagram for each transect
%                   
% Inputs:                                                                            
% - gvPtemp     = potential temperature
% - gvSalt      = salinity  
% - gvTimeTrans = time
% - AWo         = indices for AWo
% - AWr         = indices for AWr
% - LIW         = indices for LIW
% - WIW         = indices for WIW
% - WDW         = indices for WDW
% - critWIW     = criterion selected for WIW identification
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
% - figures of TS diagram with water masses in color
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 12-Apr-2018 (Melanie Juza, mjuza@socib.es)                       
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
steel=[0.27 0.51 0.71];
olive=[0.42 0.56 0.14];

% Title
str0=[critWIW,' WIW criterion - ',strrep(dataset,'_','-'),' (',datamode,')'];
str1=[section,' - ', missionName, ' - T',num2str(transectNo)];
str2=['From ' datestr(min(gvTimeTrans),1) ' to '  datestr(max(gvTimeTrans),1)];

% Figure

figure('visible',showplot); initfigall(25,20); clf;
plot(gvSalt(AWo),gvPtemp(AWo),'.','Color',burgundy); hold all 
plot(gvSalt(AWr),gvPtemp(AWr),'.','Color',steel); 
plot(gvSalt(WIW),gvPtemp(WIW),'.','Color',olive); 
plot(gvSalt(LIW),gvPtemp(LIW),'.','Color',[.5 .5 .5]); 
plot(gvSalt(WDW),gvPtemp(WDW),'.','Color','k'); hold all  
limX = [max(36,floor(min(min(gvSalt))*10)/10) min(38.8,ceil(max(max(gvSalt))*10)/10)];
limY = [max(12,floor(min(min(gvPtemp))*10)/10) ceil(max(max(gvPtemp))*10)/10];
xlim(limX);
ylim(limY);
text(limX(2)-(limX(2)-limX(1))/6,limY(2)-(limY(2)-limY(1))/20,  'AWo', 'Color',burgundy,'FontSize',font);
text(limX(2)-(limX(2)-limX(1))/6,limY(2)-(limY(2)-limY(1))/10,  'AWr', 'Color',steel,'FontSize',font);
text(limX(2)-(limX(2)-limX(1))/6,limY(2)-(limY(2)-limY(1))/6.5 ,'WIW', 'Color',olive,'FontSize',font);
text(limX(2)-(limX(2)-limX(1))/6,limY(2)-(limY(2)-limY(1))/5,   'LIW', 'Color',[.5 .5 .5],'FontSize',font);
text(limX(2)-(limX(2)-limX(1))/6,limY(2)-(limY(2)-limY(1))/4,   'WMDW','Color','k','FontSize',font);

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
imageName = ['TSdiagram_watermass_',critWIW,'_',dataset,'_',section,'_',missionName,'_T',num2str(transectNo),'_',datamode];
print('-dpng',fullfile(imageDir,imageName));
display(['    - ',imageName,'.png saved in ',imageDir])

done ='ok';

end


 
 








