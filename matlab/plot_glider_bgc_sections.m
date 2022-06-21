function [done] = plot_glider_bgc_sections(data, datamode,dataset,missionName,imageDir,logo,showplot)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function plot_glider_bgc_sections(data, datamode,dataset missionName,imageDir,logo,showplot)  
%                                                                     
% Purpose:                                                                           
% - Plot vertical sections of chlorophyll, oxygen concentration and saturation, turbidity for all the mission
%                   
% Inputs:                                                                            
% - data 
% - datamode
% - dataset
% - mission name                                       
% - figureDir = figure output directory      
     % - logo      = 'socib','copernicus','nologo'                                                                              
% - showplot  = 'on' or 'off' (to show plot on the screen)                                                                                                             
%   
% Outputs:
% - figure of vertical sections of chlorophyll with bathymetry line
% - figure of vertical sections of oxygen concentration with bathymetry line
% - figure of vertical sections of oxygen saturation with bathymetry line
% - figure of vertical sections of turbidity with bathymetry line
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: Nov-2017 (Melanie Juza, mjuza@socib.es)                       
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure
font=20;
if ~strcmp (logo, 'nologo')
im=get_logo([logo]);
end

%%% Bathymetry

% Read bathymetry
[blon blat bathy]=get_bathymetry_Smith_Sandwell;

% Interpolation of bathy on lon/lat from glider and remap in function of glider size (ntim,ndep)
bathy_interp=interp2(blon,blat,bathy,data.longitude,data.latitude);

%%% Figures of vertical sections for all the mission

for ivar=1:4

  if ivar==1
    % Chlorophyll
    datavar=data.chlorophyll;
    unit='mg.m-3';
    str0=['Chlorophyll - ',strrep(dataset,'_','-'),' (',datamode,')'];
    varname='Chl';
  elseif ivar==2
    % Oxygen concentration
    datavar=data.oxygenconc;
    unit='umol.l-1';
    str0=['Oxygen concentration - ',strrep(dataset,'_','-'),' (',datamode,')'];
    varname='O2conc';
  elseif ivar==3
    % Oxygen saturation
    datavar=data.oxygensat;
    unit='';
    str0=['Oxygen saturation - ',strrep(dataset,'_','-'),' (',datamode,')'];
    varname='O2sat';
  elseif ivar==4
    % Turbidity
    datavar=data.turbidity;
    unit='NTU';
    str0=['Turbidity - ',strrep(dataset,'_','-'),' (',datamode,')'];
    varname='Turb';
  end

  clear idx               
  idx = any(~isnan(data.ptemp),1);
  vmin=prctile(datavar(:),1);
  vmax=prctile(datavar(:),99);

  [mytime2D,mydep2D]=meshgrid(data.time(idx),data.depth);
  datavar2D=datavar(:,idx);

  if sum(idx)>100 & length(datavar(isfinite(datavar)))>0
  % Plot vertical section
  figure('visible',showplot); initfigall(40,20); clf;
  set(gca,'FontSize',font);
  %h=pcolor(data.time(idx),-data.depth, datavar(:,idx));hold on; shading flat; 
  scatter(mytime2D(:),-mydep2D(:),10,datavar2D(:),'filled');
  % Bathy contour
  hold on; plot(data.time(idx),bathy_interp(idx),'Color',[.8 .8 .8],'linewidth',2);
  fill( [data.time(idx) fliplr(data.time(idx))],  [bathy_interp(idx) fliplr(-1100*ones(1,length(bathy_interp(idx))))], [.8 .8 .8]);
  % Grid and labels
  h=colorbar;
  set(get(h, 'xLabel'), 'String',unit, 'Fontsize',font);
  set(gca,'fontsize',font);
  caxis([vmin vmax]); AA=axis; axis([AA(1) AA(2) min(bathy_interp(idx)) 0]);
  xlabel('Time','fontsize',font);
  ylabel('Depth (m)','fontsize',font);
  set(gca,'xticklabel',[])
  grid on ;
  str1=['Balearic Channels - ', missionName];
  str2=['From ' datestr(min(data.time(idx)),1) ' to '  datestr(max(data.time(idx)),1)];
  title(sprintf('%s\n%s\n%s',str0,str1,str2),'FontSize',font);
  AA=axis;axis([AA(1) AA(2) min(bathy_interp) 0]);
  % Add logo
  if ~strcmp (logo, 'nologo')
  %axes('position',[0.05,0.89,0.1,0.1])
  %imshow(im)
  axes('position',[0.07,0.89,0.05,0.1])
  image(im)  
  h = gca; h.XAxis.Visible = 'off';
  h = gca; h.YAxis.Visible = 'off';
  end
  % Print
  imageName = [ varname,'_section_',dataset,'_', missionName, '_AllMission_',datamode];
  print('-depsc2','-painters',fullfile(imageDir,imageName));
  display(['    - ',imageName,'.eps saved in ', imageDir])
  end

  clear datavar datavar2D mytime2D mydep2D

  done = i;

end

end

