function [done] = plot_glider_ts_sections(data,datamode,dataset,missionName,imageDir,logo,showplot)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function plot_glider_ts_sections(data,datamode,dataset,missionName,imageDir,logo,showplot)  
%                                                                     
% Purpose:                                                                           
% - Plot vertical sections of temperature and salinity for all the mission
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
% - figure of vertical sections of temperature with bathymetry line
% - figure of vertical sections of salinity with bathymetry line
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 19-Jun-2018 (Melanie Juza, mjuza@socib.es)                       
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


%%% Figures of temperature & salinity vertical sections for all the mission

for ivar=1:2

  if ivar==1
    % Temperature vertical section    
    datavar=data.ptemp;
    unit='{\circ}C';
    str0=['Pot. temperature - ',strrep(dataset,'_','-'),' (',datamode,')'];
    varname='T';
  elseif ivar==2
    % Salinity vertical section    
    datavar=data.salinity;
    unit='';
    str0=['Salinity - ',strrep(dataset,'_','-'),' (',datamode,')'];
    varname='S';
  end

  clear idx               
  idx = any(~isnan(datavar),1);
  vmin=prctile(datavar(:),1);
  vmax=prctile(datavar(:),99);

  % Plot vertical section
  figure('visible',showplot); initfigall(40,20); clf;
  colormap('jet')
  set(gca,'FontSize',font);
  pcolor(data.time(idx),-data.depth, datavar(:,idx));hold on;
  shading interp; 
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
  imageName = [varname,'_section_',dataset,'_', missionName, '_AllMission_',datamode];
  print('-dpng',fullfile(imageDir,imageName));
  display(['    - ',imageName,'.png saved in ', imageDir])

  clear datavar

  done = i;

end

end

