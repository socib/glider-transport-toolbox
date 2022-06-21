function transects_vector=get_glider_transects_BC(lon,lat,section,missionName,datamode,imageDir,logo,showplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function detect_glider_transects_Balearic_Channels(lon,lat,section,missionName,datamode,imageDir,logo,showplot)
%
% Purpose:
% - Separate glider transects in a selected Balearic Channel
%
% Inputs: 
% - lon         = longitude 
% - lat         = latitude
% - section     = 'IbizaChannel' or 'MallorcaChannel'
% - missionName = mission name
% - datamode    = 'rt' o 'dt'                                                           
% - imageDir    = figure output directory  
% - logo       = 'socib','copernicus','nologo'
% - showplot    = 'on' or 'off' (to show plot on the screen)  
%                                                                                                           
% Outputs: 
% - vector containing transect numbers in Ibiza Channel (0 outside)
% - plot transect detection vs longitude
%
% Date of creation: June-2018 (adapted from get_glider_transects_IC)
%
% (Melanie Juza, Baptiste Mourre, SOCIB)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure
font=17;
if ~strcmp (logo, 'nologo')
im=get_logo([logo]);
end

strtitle=['Detected transects in ',section,' for ',missionName,' (',datamode,')'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Define Channel latitude & longitud bands
if strcmp (section, 'IbizaChannel')
  latmin=38.9; latmax=39.1;
  lonmin=0.1;  lonmax=1.13;
elseif strcmp (section, 'MallorcaChannel')
  latmin=39;   latmax=39.6;
  lonmin=1.6;  lonmax=2.3;
end

%%% Define threshold in longitude
if strcmp (section, 'IbizaChannel')
    lon_thres_east=0.95;
    lon_thres_West=0.2;
elseif strcmp (section, 'MallorcaChannel')
    lon_thres_east=2.2;
    lon_thres_West=1.7;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize vector
transects_vector=zeros(length(lon),1);

% Select indices in the Channel box
ii_BC=find(lon>=lonmin & lon<=lonmax & lat>=latmin & lat<=latmax);

if ~isempty(ii_BC)
    % Find "East" points with longitude > lon_thres_east
    ii_East=find(lon(ii_BC)>lon_thres_east);
    % Detect jumps, corresponding to indices when the glider crosses lon_thres_east moving towards the West (adding first and last incices)
    ii_East_jump=[1;find(diff(ii_East)>1);length(ii_East)];
    % Consider mean index when the glider is East of lon_thres_east, and store it in ii_East_transects (indices of end/start of transects in East side)
    ii_East_transects=[];
    for ki=1:length(ii_East_jump)-1
        ii_East_transects(ki)=floor(mean(ii_East(ii_East_jump(ki)+1:ii_East_jump(ki+1))));
    end
    if ~isnan(ii_East_transects)>0
        % Do the same with the western side     
        ii_West=find(lon(ii_BC)<lon_thres_West);
        % Detect jumps, corresponding to indices when the glider crosses lon_thres_West moving towards East (adding first and last incices)
        ii_West_jump=[1;find(diff(ii_West)>1);length(ii_West)];
        % Consider mean index when the glider is West of lon_thres_West, store it in ii_West_transects (indices of end/start of transects in West side)
        ii_West_transects=[];
        for ki=1:length(ii_West_jump)-1
            ii_West_transects(ki)=floor(mean(ii_West(ii_West_jump(ki)+1:ii_West_jump(ki+1))));
        end
        if ~isnan(ii_West_transects)>0
            % Combine all end/start of transects
            ii_change_transects=sort([ii_East_transects ii_West_transects]);
            % Number transects, keep 0 outside the Channel
            kkt=0;
            for kt=1:length(ii_change_transects)-1
                if ~isempty(find(lon(ii_BC(ii_change_transects(kt)):ii_BC(ii_change_transects(kt+1)))<=lon_thres_West)) ...
                   &  ~isempty(find(lon(ii_BC(ii_change_transects(kt)):ii_BC(ii_change_transects(kt+1)))>=lon_thres_east))
                    kkt=kkt+1;
                    transects_vector(ii_BC(ii_change_transects(kt)):ii_BC(ii_change_transects(kt+1)))=kkt;
                end
            end          
            disp([num2str(kt) ' transect(s) detected in ' section])
        else
            disp(' No data close to the western Channel')
            ii_change_transects=[];
        end
    else
        disp(' No data close to the eastern Channel')
        ii_change_transects=[];
    end
else
    disp([' No data in ',section])
    ii_change_transects=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot full longitude evolution and selected Balearic Channel transects changes    
nt=length(unique(transects_vector))-1;
figure('visible',showplot); initfigall(30,20); clf;
set(gca,'FontSize',font)
plot(lon,'LineWidth',2);hold on
for kt=1:nt
  i1=min(find(transects_vector==kt));
  i2=max(find(transects_vector==kt));
  plot([i1 i1],ylim,'r')
  hold on
  if (kt==1)
    plot([i1 i1],ylim,'r','LineWidth',3)
  end
  if (kt==nt)
    plot([i2 i2],ylim,'r','LineWidth',3)
  end
end
xlabel('Index in glider data','FontSize',font)
ylabel('Longitude','FontSize',font)
grid on;
title (strtitle,'fontsize',font)
% Logo
if ~strcmp (logo, 'nologo')
  if strcmp (logo, 'copernicus')
    axes('position',[0.01,0.92,0.09,0.09])
  elseif strcmp (logo, 'socib')
    %axes('position',[0.01,0.89,0.1,0.1])
    axes('position',[0.025,0.89,0.065,0.1])
  end
  image(im)  
  h = gca; h.XAxis.Visible = 'off';
  h = gca; h.YAxis.Visible = 'off';
end
% Print
imageName = ['DetectedTransects_',missionName,'_',datamode,'_',section];
print('-dpng',fullfile(imageDir,imageName)); 
display(['    - ',imageName,'.png saved in ', imageDir])


 
end

