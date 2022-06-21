function [datafilename]=processing_transect_glider(filename,datamode,datatype,dataset,section,matfileDir,figureDir,optionbgc,optionfile,optionplot,logo,showplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function processing_transect_glider(filename,datatype,section,matfileDir,     
%                           figureDir,option.bgc,optionsavefile,optionplot,showplot) 
%                                                                                    
% Inputs:                                                                            
% - filename  = glider file (L1 product)              
% - datamode  ='rt' o 'dt'                                                           
% - datatype  = 'G' (glider)                                
% - dataset   = 'glider'                      
% - section   = 'IbizaChannel'                                                       
% - matfileDir= data output directory                                                
% - figureDir = figure output directory                                              
% - optionbgc=  true or false (to process biogeochemical data)                       
% - optionfile= true or false (to save file)                                          
% - optionplot= true or false (to create figures)  
% - logo      = 'socib','copernicus','nologo'                                  
% - showplot  = 'on' or 'off' (to show plot on the screen)                           
%                                                                                    
% General purpose:                                                                   
% - Quality control processing, vertical interpolation, transect detection           
%                                                                                    
% Purpose:                                                                           
% - Discards profiles with large gaps if numel(goodRows)>dataMin &                   
%                                     if max(abs(diff(CTD_data(goodRows, 2))))<dbGap 
% - Interpolate profiles in vertical & takes mean loc, mean time as representitative 
% - Deletes known erroneous profiles                                                 
% - Detects the transects along the channel                                          
%                                                                                    
% Outputs:                                                                           
% - Save TS data in structure                                                       
% - Plots: TS diagram, TS (& BGC) vertical sections, transect maps (optional)                            
%                                                                                    
%                                                                                    
% Adapted from SOCIB glider toolbox (Emma Heslop, Glider Facility, Data Centre, 2015)
%                                                                                                     
% Author: Melanie Juza - SOCIB                                                       
%         mjuza@socib.es                                                             
%                                                                                    
% Last modification: 20-Jun-2022)                                         
%                                                                                    
% (History: modifed version v8, no timeseries/spectra check - glider only)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off all

disp(['1. Data processing',]);
disp(['   Start time: ', datestr((now), 'yyyy-mm-dd HH:MM:SS')]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set QC 
dataMin   = 10;  % profile must have 10 data points
dbGap     = 5;   % profile must have data gap < 5 db
minLength = 20 ; % profile must length of data > 20 db

% Depth resolution (in m) for interpolation
depthResolution = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Variables
  
% Read variables
fileContent = nc_info(filename);
for varIdx = 1:length(fileContent.Dataset)
  varName = fileContent.Dataset(varIdx).Name;                    
  eval([varName, ' = nc_varget(filename, ''', varName, ''');']); 
  eval([varName,'(',varName,'>1e+36)=NaN;']);
end;
clear fileContent varIdx varName;
  
% Time (automatically select ctd time if exists)
if exist('time_ctd', 'var') 
  timeElected = time_ctd;
else
  timeElected = time;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output file and directory

timeElected(timeElected==0)=NaN;
qctime=timeElected(find(isfinite(timeElected)));
start=qctime(1)/(60*60*24) + datenum(1970,1,1,0,0,0);
missionName=['canales',datestr(start,'mmmyyyy')];

% Make folders for images
matoutDir= fullfile(matfileDir,missionName);
mkdir(matoutDir);
imageDir = fullfile(figureDir,missionName);
mkdir(imageDir);
imageDir2= fullfile(imageDir,section);
mkdir(imageDir2);

%%% Output filename
if optionbgc
  ext='_bgc.mat';
else
  ext='.mat';
end
datafilename = fullfile(matoutDir,[dataset,'_',missionName,'_transportGV_',datamode,'_',section,ext]);

if exist(datafilename)
  display(['   ',datafilename,' already exists'])
  load(datafilename);
  if exist('dataTS')  
    display('   --> dataTS have been already reprocessed for transport computation')
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization

% Number of profiles
C = unique(profile_index);
nbProfiles = length(C);

% Max and Min depths for bin and interpolation
maxDepth    = round(max(pressure));
minDepth    = 0; 
depthLimits = [minDepth, maxDepth];
depthRange  = (minDepth:depthResolution:maxDepth);
depthLevels = length(depthRange);

% Create matrixes
posMatrix      = nan(4,          nbProfiles); % 4 rows: lon, lat, 2 x time (seconds from 1970-01-01,00:00:00)
maxDepthMatrix = nan(1,          nbProfiles); % max depth
minDepthMatrix = nan(1,          nbProfiles); % min depth
tempMatrix     = nan(depthLevels,nbProfiles); % temp matrix
saltMatrix     = nan(depthLevels,nbProfiles); % salinity matrix 
if optionbgc
  chloMatrix   = nan(depthLevels,nbProfiles); % chlorophyll matrix
  oconMatrix   = nan(depthLevels,nbProfiles); % oxygen concentration matrix 
  osatMatrix   = nan(depthLevels,nbProfiles); % oxygen saturation matrix 
  turbMatrix   = nan(depthLevels,nbProfiles); % turbidity matrix 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% data processing for each profile

erroneusProfiles = [];
  
for prfIdx = 1:nbProfiles  

  clear CTD_data profilePres origProfileTemp profileSalt    
  indexRange = find(profile_index == C(prfIdx));

  % Variables
  if strcmp ( datatype, 'G') 
      CTD_data = [timeElected(indexRange), pressure(indexRange), temperature(indexRange), salinity(indexRange)];
      if optionbgc
        CTD_data = [timeElected(indexRange), pressure(indexRange), temperature(indexRange), salinity(indexRange), ...
                    chlorophyll(indexRange), oxygen_concentration(indexRange), oxygen_saturation(indexRange), turbidity(indexRange)];
      end
  end

  % Find rows where 4 variables available (time, pressure, temperature, salinity)
  goodRows = find(sum(isnan(CTD_data(:,1:4)), 2) == 0 );

  %% Discard profiles if number of data<dataMin & gaps>dbGap
  if numel(goodRows) > dataMin && max(abs(diff(CTD_data(goodRows, 2)))) < dbGap && abs(CTD_data(goodRows(1), 2)-CTD_data(goodRows(end), 2)) > minLength

    CTD_data    = CTD_data(goodRows,:);
    profilePres = CTD_data(:,2);
    profileTemp = CTD_data(:,3); 
    profileSalt = CTD_data(:,4);
    if optionbgc
      profileChlo = CTD_data(:,5);
      profileOcon = CTD_data(:,6); 
      profileOsat = CTD_data(:,7);
      profileTurb = CTD_data(:,8);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Interpolate in vertical (bin and interpolate to regular depth in vertical)

    % Normally 1 m resolution (approx. glider sampling)
    clear interpRange interpprofileTemp interpprofileSalt
    minProfDepth      = floor(min(profilePres));
    maxProfDepth      = ceil(max(profilePres));
    interpRange       = find (depthRange >=minProfDepth & depthRange <= maxProfDepth); 
    interpprofileTemp = NaN(1, length(depthRange));
    interpprofileSalt = NaN(1, length(depthRange));
    
    % Remove duplicated pressure in profile before interpolation
    p = find (diff(profilePres) == 0);
    if ~isempty(p) 
      n              = (1: length(profilePres));
      noDoublesIndex = setdiff(n,p);
      profilePres    = profilePres(noDoublesIndex);
      profileTemp    = profileTemp(noDoublesIndex);
      profileSalt    = profileSalt(noDoublesIndex);
    end 
    
    % Interpolate in Vertical 
    interpprofileTemp(interpRange) = interp1(profilePres, profileTemp, depthRange(interpRange), 'linear', NaN)';
    interpprofileSalt(interpRange) = interp1(profilePres, profileSalt, depthRange(interpRange), 'linear', NaN)';
      
    % Introduce profile information into a matrix of the transect
    tempMatrix(:,prfIdx) = interpprofileTemp(:); 
    saltMatrix(:,prfIdx) = interpprofileSalt(:);
    
    % Assume profile is vertical at mean of all lat, lon, time in profile  
    posMatrix(1,prfIdx) = nanmean(latitude(indexRange));
    posMatrix(2,prfIdx) = nanmean(longitude(indexRange));
    posMatrix(3,prfIdx) = nanmean(timeElected(indexRange));
    
    profileMatrix(prfIdx) = C(prfIdx);
   
    % Create max/min depth matrix 
    maxDepthMatrix(prfIdx) = maxProfDepth;
    minDepthMatrix(prfIdx) = minProfDepth;

    if optionbgc
      clear interpprofileChlo interpprofileOcon interpprofileOsat interpprofileTurb
      interpprofileChlo = NaN(1, length(depthRange));
      interpprofileOcon = NaN(1, length(depthRange));
      interpprofileOsat = NaN(1, length(depthRange));
      interpprofileTurb = NaN(1, length(depthRange));    
      if ~isempty(p) 
        profileChlo = profileChlo(noDoublesIndex);
        profileOcon = profileOcon(noDoublesIndex);
        profileOsat = profileOsat(noDoublesIndex);
        profileTurb = profileTurb(noDoublesIndex);
      end
      if length(profilePres)==length(profileChlo)
        interpprofileChlo(interpRange) = interp1(profilePres, profileChlo, depthRange(interpRange), 'linear', NaN)';
        interpprofileOcon(interpRange) = interp1(profilePres, profileOcon, depthRange(interpRange), 'linear', NaN)';
        interpprofileOsat(interpRange) = interp1(profilePres, profileOsat, depthRange(interpRange), 'linear', NaN)';
        interpprofileTurb(interpRange) = interp1(profilePres, profileTurb, depthRange(interpRange), 'linear', NaN)';
      end
      chloMatrix(:,prfIdx) = interpprofileChlo(:); 
      oconMatrix(:,prfIdx) = interpprofileOcon(:);
      osatMatrix(:,prfIdx) = interpprofileOsat(:); 
      turbMatrix(:,prfIdx) = interpprofileTurb(:);
    end
   
  else
    erroneusProfiles = [erroneusProfiles, prfIdx];
  end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Processing interpolated profiles

% Remove erroneous profiles from the matrices
if length(erroneusProfiles)>=1
  posMatrix      = posMatrix(:,     setdiff(1:nbProfiles, erroneusProfiles));
  maxDepthMatrix = maxDepthMatrix(:,setdiff(1:nbProfiles, erroneusProfiles));
  minDepthMatrix = minDepthMatrix(:,setdiff(1:nbProfiles, erroneusProfiles));
  tempMatrix     = tempMatrix(:,    setdiff(1:nbProfiles, erroneusProfiles));
  saltMatrix     = saltMatrix(:,    setdiff(1:nbProfiles, erroneusProfiles));
  if optionbgc
    chloMatrix   = chloMatrix(:,    setdiff(1:nbProfiles, erroneusProfiles));
    oconMatrix   = oconMatrix(:,    setdiff(1:nbProfiles, erroneusProfiles));
    osatMatrix   = osatMatrix(:,    setdiff(1:nbProfiles, erroneusProfiles));
    turbMatrix   = turbMatrix(:,    setdiff(1:nbProfiles, erroneusProfiles));
  end
  profileMatrix  = profileMatrix(:, setdiff(1:nbProfiles, erroneusProfiles));
  disp(['   Profiles in mission (good/total): ', num2str(length(C)-length(erroneusProfiles)),'/',num2str(length(C))]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save data in structure 

dataTS.sciTime     = posMatrix(3, :);
dataTS.latitude    = posMatrix(1, :);
dataTS.longitude   = posMatrix(2, :);
dataTS.profile     = profileMatrix;
dataTS.temperature = tempMatrix;
dataTS.salinity    = saltMatrix;
dataTS.depth       = depthRange;
dataTS.maxDepth    = maxDepthMatrix;
dataTS.minDepth    = minDepthMatrix; 
if optionbgc
  dataTS.chlorophyll = chloMatrix;
  dataTS.oxygenconc  = oconMatrix;
  dataTS.oxygensat   = osatMatrix;
  dataTS.turbidity   = turbMatrix;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Specific corrections

% Remove possible temperature spikes in glider data
if strcmp (datatype, 'G')
  tSpikes = find(dataTS.temperature < 10);
  dataTS.temperature(tSpikes) = NaN;
  dataTS.salinity(tSpikes)    = NaN;
  sSpikes = find(dataTS.salinity < 35);
  dataTS.temperature(sSpikes) = NaN;
  dataTS.salinity(sSpikes)    = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Additional variables

% Additional variables (pot. temperature, density) 
pressMatrix  = dataTS.depth(:) * ones(1,size(dataTS.temperature,2)); 
if strcmp ( datatype, 'G') 
  dataTS.ptemp   = sw_ptmp(dataTS.salinity, dataTS.temperature, pressMatrix, 0); 
elseif strcmp ( datatype, 'M') 
  dataTS.ptemp   = dataTS.temperature; 
end
dataTS.density = sw_dens(dataTS.salinity, dataTS.ptemp, 0);
dataTS.time    = dataTS.sciTime/(60*60*24) + datenum(1970,1,1,0,0,0) ; % sciTime to matlab time
clear tempMatrix saltMatrix posMatrix profileMatrix maxDepthMatrix minDepthMatrix  

% Glider characteristics
dataTS.mission =missionName;
if datatype == 'G'
  words = regexp(strrep(filename,'/',' '),'\s+', 'split');
  file=words{end};
  dataTS.gliderType=file(9:15);
elseif datatype == 'M'
  dataTS.gliderType='model';
end

% Get transects in the Channel 
transects_vector=get_glider_transects_BC(dataTS.longitude',dataTS.latitude',section,missionName,datamode,imageDir2,logo,showplot);
dataTS.transect = transects_vector;   
transects_vector_BC=transects_vector(transects_vector~=0);
if ~isempty(transects_vector_BC)
  indtr=[1;find(diff(transects_vector_BC))];
  if strcmp (section, 'IbizaChannel')
    dataTS.ibiza = transects_vector_BC(indtr+1); 
  elseif strcmp (section, 'MallorcaChannel')
    dataTS.mallorca = transects_vector_BC(indtr+1); 
  end
else
  if strcmp (section, 'IbizaChannel')
    dataTS.ibiza = NaN;
  elseif strcmp (section, 'MallorcaChannel')
    dataTS.mallorca = NaN; 
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots

%%% Save file
if optionfile
  save(datafilename, 'dataTS')
  disp(['    - ',dataset,'_',missionName,'_transportGV_',datamode,'_',section,ext ' saved in ', matoutDir]);
end

%%% Plots all transects of the mission

% Temperature/salinity diagram
if optionplot
  [done] = plot_glider_ts_diagram_IC_MC(dataTS,datamode,dataset,missionName,imageDir,logo,showplot);
end

% Temperature and salinity vertical sections
if optionplot
  [done] = plot_glider_ts_sections(dataTS,datamode,dataset,missionName,imageDir,logo,showplot);
  if optionbgc
  [done] = plot_glider_bgc_sections(dataTS,datamode,dataset,missionName,imageDir,logo,showplot);
  end
end

% Map all transects  
if optionplot
  [done] = plot_glider_trajectory(dataTS.longitude,dataTS.latitude,missionName,datamode,imageDir,logo,showplot);
end

% Stop deployment processing logging.
disp(['   End time: ' , datestr((now), 'yyyy-mm-dd HH:MM:SS')]);




