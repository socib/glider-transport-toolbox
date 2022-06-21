function [datafilename]=processing_transport_glider(datafilename,datamode,datatype,dataset,section,smoothParam,critWIW,matfileDir,figureDir,optionbgc,optionfile,optionplot,logo,showplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function processing_transport_glider(datafilename,datamode,datatype,section, 
%      smoothParam,matfileDir,figureDir,optionbgc,optionfile,optionplot,logo,showplot )  
%                                                                                    
% Inputs:                                                                            
% - datafilename = processed data for glider (from processing_transect_glider.m)                   
% - datamode     = 'rt' o 'dt'                                                       
% - datatype     = 'G' (glider)                                  
% - dataset      = 'glider'                         
% - section      = 'IbizaChannel'                                                    
% - smoothParam  = smoothing parameter (in km)                                       
% - critWIW      ='fixed-range','visual','geometric'                                 
% - matfileDir   = data output directory                                             
% - figureDir    = figure output directory                                           
% - optionbgc    = true or false (to process biogeochemical data)                       
% - optionfile   = true or false (to save file)                                       
% - optionplot   = true or false (to create figures) 
% - logo         = 'socib','copernicus','nologo'                         
% - showplot     = 'on' or 'off' (to show plot on the screen)                        
%                                                                                    
% General purpose:                                                                   
% - Data processing for transport computation                                        
%                                                                                    
% Purpose:                                                                           
% - Projection to standard line                                                      
% - Fill surface values                                                              
% - Bin and interpolate in vertical                                                  
% - Bin and interpolate in horizontal                                                
% - Fill to channel (defined with bathymetry)                                        
% - Smooth in horizontal (moving average)                                            
% - Compute GV calc per profile pair (common depth is zero ref)                                              
% - Compute transports: total N & S, and by water mass (for channel)                 
%                                                                                    
% Outputs:                                                                           
% - Saves GV and transports data in structure                                        
% - Plots: GV vertical sections (optional)                                           
%                                                                                    
%                                                                                    
% Adapted from SOCIB glider toolbox (Emma Heslop, Glider Facility, Data Centre, 2015)
%                                                                                                     
% Author: Melanie Juza - SOCIB                                                       
%         mjuza@socib.es                                                             
%                                                                                    
% Last modification: 20-Jun-2022                          
%                                                                                    
% History: modified version v11 (only gliders)                 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off all
 
disp(['2. Transport processing']);
disp(['   Start time: ', datestr((now), 'yyyy-mm-dd HH:MM:SS ')]);

if exist(datafilename)
  load(datafilename);
  if exist('dataGV') & exist('transports') & exist('statistics')
    display('   --> dataGV & transports & statistics have been already computed')
    return
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Options for bin/interp
options.binVert    = true;  % Bin in vertical before interpolation
binResolutionVert  = 5;     % Vertical resolution of T/S data bins (m) 
binResolutionHoriz = 2;     % Horizontal resolution of T/S data bins km (km)

% Define lat and lon of standard transect line (glider way points)
if ismember('IbizaChannel', section) 
  latPts=[38.985 38.985];  
  lonPts=[ 0.06   1.16 ]; 
elseif ismember('MallorcaChannel', section) 
  latPts=[39.222 39.5007]; 
  lonPts=[01.653 02.182 ];   
end  

% Options for gv calc
options.gvGridOrig = true;  % gv interp back to orig. horiz. bin/interp positions or keep mean positions from gv calc.
options.useRefVel  = false; % Use a reference level velocity for gv calculation (false: ref velocity = 0)
useRefVels         = 1;

% Options channel 
options.channelTransports = true; % To save only channel transports in transport file, if not whole transect (channel and shelf)
shelfDepth                = 200;  % Minimum depth of profile data for main channel  

% Parameters
limit       = 15;         % Limit alert for profile distance to line (km)
dbar2Pascal = 10000;      % Factor conversion: from decibars to Pascals for dynamic Ht
omega       = 7.29*10^-5; % s^-1
g           = 9.81;       % m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read file

  % Data
  load(datafilename);
  missionName=dataTS.mission;

  % Transects  
  if strcmp(section, 'IbizaChannel') 
    transectList = dataTS.ibiza;
  elseif strcmp(section, 'MallorcaChannel') 
    transectList = dataTS.mallorca;
  end
  disp(['   Transect numbers in ',section,': ', num2str(transectList')]);

  % Dynamic height reference level for transport calculation
  clear depthIdx depth depthDynHeight
  minDynHeightRefLevel = 0; 
  maxDynHeightRefLevel = max(dataTS.depth); 
  depthRange = dataTS.depth;
  depthIdx = find(depthRange >= minDynHeightRefLevel & depthRange <= maxDynHeightRefLevel); 
  depthDynHeight = depthRange(depthIdx);                                                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialization

  % Set up files
  gvLat   = [];   gvTransport = [];
  gvLon   = [];   gvDist      = [];
  gvTime  = [];   gvTransect  = [];
  gvDepth = [];   tsPtemp     = [];
  gvPtemp = [];   tsSalt      = [];
  gvSalt  = [];   
  gvDens  = [];   
  gv      = [];
  if optionbgc
   gvChlor = [];
   gvOconc = [];
   gvOsatu = [];
   gvTurbi = [];
  end
  gvBathy  = [];

  % Transport calculation & plotting program 
  count      = 1; % counts for plots
  countFinal = 1; % counts for final plots   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output directory

  matoutDir= fullfile(matfileDir,missionName);
  imageDir = fullfile(figureDir,missionName,section);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Transect processing...
  
  if length(transectList)>=1 & ~isnan(transectList)

    for t=1:length(transectList) 

      transectNo = transectList(t);
      clear transectInd 
      transectInd = find(dataTS.transect==transectNo); 
      disp(['   T',num2str(transectNo),': number profiles is ' , num2str(length(transectInd))]);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% Projection of lon/lat on standard line

      % Orientate all data in W to E direction
      if  dataTS.longitude(1,transectInd(end)) < dataTS.longitude(1,transectInd(1)) 
        transectInd = flipdim(transectInd,2);
      end

      % Project lat/lon onto standard lines 
      clear latProj lonProj dist2line dist2Prof standardLine 
      [lonProj, latProj, dist2Prof, gridLine] = stdLineProjection(dataTS, transectInd, latPts, lonPts, binResolutionHoriz);  
      disp(['   Transect length: ', num2str(max(dist2Prof), 2), 'km'])
      % Distance to standard line 
      dist2Line = distance([dataTS.latitude(transectInd)';latProj(:,1)],[dataTS.longitude(transectInd)';lonProj(:,1)],'km');

      % Remove duplication of profiles
      [b, p, n] = unique(dist2Prof);
      if length(b)<length(transectInd)
        transectInd = transectInd(p);
        dist2Prof   = dist2Prof(p);
        lonProj     = lonProj(p);
        latProj     = latProj(p);
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% Fill gaps

      % Fill surface different methods for different platforms
      clear fillProfTemp fillProfSalt fillProfChlo fillProfOcon fillProfOsat fillProfTurb ptemp salinity chlorophyll oxygenconc oxygensat turbidity
      fillProfTemp = dataTS.ptemp(depthIdx,transectInd);
      fillProfSalt = dataTS.salinity(depthIdx,transectInd);
      if optionbgc
        fillProfChlo = dataTS.chlorophyll(depthIdx,transectInd);
        fillProfOcon = dataTS.oxygenconc(depthIdx,transectInd);
        fillProfOsat = dataTS.oxygensat(depthIdx,transectInd);
        fillProfTurb = dataTS.turbidity(depthIdx,transectInd);
      end

      inst=dataTS.gliderType;
      [gliderType] = findGliderType(inst);

      if strcmp(gliderType, 'g1') || strcmp(gliderType, 'g2') || strcmp(gliderType, 'g3')
        
        % Fill surface with horizontal (linear) interpolation between all values & profiles that go to the surface
        clear surfaceFillDepth
        surfaceFillDepth = max(dataTS.minDepth) +1;
        surfaceFillDepth = find(dataTS.depth == surfaceFillDepth);  
        for ii=1:surfaceFillDepth
          clear realCols
          realCols = find(~isnan(fillProfTemp(ii,:)));
          if length(realCols)>=2
            fillProfTemp(ii,:) = interp1(dist2Prof(realCols),fillProfTemp(ii,realCols),dist2Prof,'pchip', NaN);   
            fillProfSalt(ii,:) = interp1(dist2Prof(realCols),fillProfSalt(ii,realCols),dist2Prof,'pchip', NaN);   
          end
        end            
        % Fill any parts in surface that outside the interpolation      
        for ii=1:surfaceFillDepth 
          clear realCols
          realCols = find(~isnan(fillProfTemp(ii,:)));
          if length(realCols)>=2 
            fillProfTemp(ii,:) = interp1(dist2Prof(realCols),fillProfTemp(ii,realCols),dist2Prof,'nearest','extrap');   
            fillProfSalt(ii,:) = interp1(dist2Prof(realCols),fillProfSalt(ii,realCols),dist2Prof,'nearest','extrap');      
          end
        end
        % Fill matrix 0m depth
        surfaceFillDepth = 0;
        surfaceFillDepth = find(dataTS.depth == surfaceFillDepth);
        fillProfTemp(1:surfaceFillDepth,:) = fillProfTemp(surfaceFillDepth+1,:);

      elseif ismember('sg', gliderType) 

        % Delete surface 5 m and fill with last (highest value) 
        surfaceFillDepth = 5;
        surfaceFillDepth = find(dataTS.depth == surfaceFillDepth);
        for ii=1:surfaceFillDepth
          fillProfTemp(ii,:) = fillProfTemp(surfaceFillDepth+1,:);
          fillProfSalt(ii,:) = fillProfSalt(surfaceFillDepth+1,:);
        end

      end

      % Fill any other gaps, Slocum where glider surfaced (nearest neighbour)
      for ii=1:length(depthDynHeight)
        realCols = find(~isnan(fillProfTemp(ii,:)));
        if length(realCols)>=2
          fillProfTemp(ii,:) = interp1(dist2Prof(realCols),fillProfTemp(ii,realCols),dist2Prof,'nearest', NaN);   
          fillProfSalt(ii,:) = interp1(dist2Prof(realCols),fillProfSalt(ii,realCols),dist2Prof,'nearest', NaN);      
        end
      end       

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% Bins and interpolations

      % Bin data in vertical and interpolate
      if options.binVert
        [binProfTemp,binProfSalt,binDepths] = binVert(fillProfTemp,fillProfSalt,dataTS.depth,binResolutionVert,max(dataTS.depth),min(dataTS.depth));
        [interpProfTemp1,interpProfSalt1,interpDepths] = interpVert(binProfTemp,binProfSalt,binDepths,binResolutionVert,max(dataTS.depth),min(dataTS.depth),'linear',NaN);
      else  
        [interpProfTemp1,interpProfSalt1,interpDepths] = interpVert(fillProfTemp,fillProfSalt,dataTS.depth,binResolutionVert,max(dataTS.depth),min(dataTS.depth),'linear',NaN);
      end
      depthLevels = interpDepths;
      interpProfTemp1(1,:) = interpProfTemp1(2, :);
      interpProfSalt1(1,:) = interpProfSalt1(2, :);
      if optionbgc
        if options.binVert
          [binProfChlo,binProfOcon,binProfOsat,binProfTurb,binDepths] = binVert_bgc(fillProfChlo,fillProfOcon,fillProfOsat,fillProfTurb,dataTS.depth,binResolutionVert,max(dataTS.depth),min(dataTS.depth));
          [interpProfChlo1,interpProfOcon1,interpProfOsat1,interpProfTurb1,interpDepths] = interpVert_bgc(binProfChlo,binProfOcon,binProfOsat,binProfTurb,binDepths,binResolutionVert, ...
                                                                                                          max(dataTS.depth),min(dataTS.depth),'linear',NaN);
        else
          [interpProfChlo1,interpProfOcon1,interpProfOsat1,interpProfTurb1,interpDepths] = interpVert_bgc(fillProfChlo,fillProfOcon,fillProfOsat,fillProfTurb,dataTS.depth,binResolutionVert,...
                                                                                                          max(dataTS.depth),min(dataTS.depth),'linear',NaN);
        end
        depthLevels = interpDepths;
        interpProfChlo1(1,:) = interpProfChlo1(2, :);
        interpProfOcon1(1,:) = interpProfOcon1(2, :);
        interpProfOsat1(1,:) = interpProfOsat1(2, :);
        interpProfTurb1(1,:) = interpProfTurb1(2, :);
      end

      % Bin data in horizontal and interpolate
      clear binProfTemp binProfSalt interpProfTemp2 interpProfSalt2 
      clear binDist 
      [binProfTemp,binProfSalt,binProfileTime,binDist] = binHoriz(interpProfTemp1,interpProfSalt1,dataTS.time(transectInd),dist2Prof,depthLevels,gridLine,binResolutionHoriz);
      [interpProfTemp2,interpProfSalt2] = interpHoriz(binProfTemp,binProfSalt,binDist,'linear',NaN);
      if optionbgc
        clear binProfChlo binProfOcon binProfOsat binProfTurb interpProfChlo2 interpProfOcon2 interpProfOsat2 interpProfTurb2 
        clear binDist 
        [binProfChlo,binProfOcon,binProfOsat,binProfTurb,binProfileTime,binDist] = binHoriz_bgc(interpProfChlo1,interpProfOcon1,interpProfOsat1,interpProfTurb1, ...
                                                                                                dataTS.time(transectInd),dist2Prof,depthLevels, gridLine, binResolutionHoriz);
        [interpProfChlo2,interpProfOcon2,interpProfOsat2,interpProfTurb2] = interpHoriz_bgc(binProfChlo,binProfOcon,binProfOsat,binProfTurb,binDist,'linear',NaN);
      end

      % Find lat and lon for bin point locations
      clear latBin lonBin
      az = distance([latPts(1,1);latPts(1,2)],[lonPts(1,1);lonPts(1,2)],'km');
      [latBin, lonBin] = reckon(latPts(1,1),lonPts(1,1),km2deg(binDist), az);  

      % Interp time onto bin/interp grid
      clear interpProfTime
      interpProfTime(1,:) = interp1(dist2Prof, dataTS.time(transectInd), binDist, 'linear', NaN);
      clear ptemp salinity chlorophyll oxygen_con oxygen_sat turbidity standardLine
      ptemp        = interpProfTemp2;
      salinity     = interpProfSalt2;   
      if optionbgc
        chlorophyll  = interpProfChlo2;
        oxygen_con   = interpProfOcon2;
        oxygen_sat   = interpProfOsat2;
        turbidity    = interpProfTurb2;
      end
      standardLine = binDist;   

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% Fill to channel bottom before smoothing

      % Create bathymetry and interpolate to stdLine
      clear bathy segmentLon segmentLat segmentDepth interpBathy dist2ProfBathy      
      [segmentLon, segmentLat, segmentDepth] = get_bathymetry_sections(latPts , lonPts, datatype, section);
      bathy = -segmentDepth;
        %%% Deals with possible errors in bathymetry field 
        indjump=find(abs(diff(bathy))>150);
        if ~isempty(indjump)
          for i=1:length(indjump)
            bathy(indjump(i) + 1) = 0.5*(bathy(indjump(i)) + bathy(indjump(i) + 2));
          end
        end
        %%%
      dist2ProfBathy = NaN(1,length(segmentLat));
      for j=1:length(segmentLat)  
        dist2ProfBathy(j) = distance([latPts(1,1);segmentLat(1,j)],[lonPts(1,1);segmentLon(1,j)],'km');
      end
    
      % Interp bathy to match profiles along 'standard' line with bin/interp
      clear indReal interpBathy fillProfTemp2 fillProfSalt2
      indReal = find(~isnan(bathy));
      interpBathy(:,:) = interp1(dist2ProfBathy(:,indReal), bathy(:,indReal), standardLine(1,:)); 
      indReal = find(~isnan(interpBathy));
      interpBathy(:,:) = interp1( standardLine(:,indReal), interpBathy(:,indReal), standardLine(1,:),'nearest','extrap'); 

      % Fill T/S to channel bathymetry 
      fillProfTemp2 = nan(size(ptemp,1), size(ptemp,2));
      fillProfSalt2 = fillProfTemp2;
      for ii=1:length(depthLevels)
        realCols = find(~isnan(ptemp(ii,:)));   
        if length(realCols)>=2
          fillProfTemp2(ii,:) = interp1(standardLine(realCols),ptemp(ii,realCols),   standardLine, 'nearest','extrap');   
          fillProfSalt2(ii,:) = interp1(standardLine(realCols),salinity(ii,realCols),standardLine, 'nearest','extrap');          
        elseif length(realCols)==1
          fillProfTemp2(ii,:) = ptemp(ii,realCols);   
          fillProfSalt2(ii,:) = salinity(ii,realCols);            
        end
      end
      if optionbgc
        fillProfChlo2 = nan(size(chlorophyll,1), size(chlorophyll,2));
        fillProfOcon2 = fillProfChlo2;
        fillProfOsat2 = fillProfChlo2;
        fillProfTurb2 = fillProfChlo2;
        for ii=1:length(depthLevels)
          realCols = find(~isnan(oxygen_con(ii,:)));   
          if length(realCols)>=2
            fillProfChlo2(ii,:) = interp1(standardLine(realCols),chlorophyll(ii,realCols),standardLine, 'nearest','extrap');   
            fillProfOcon2(ii,:) = interp1(standardLine(realCols),oxygen_con(ii,realCols), standardLine, 'nearest','extrap'); 
            fillProfOsat2(ii,:) = interp1(standardLine(realCols),oxygen_sat(ii,realCols), standardLine, 'nearest','extrap');   
            fillProfTurb2(ii,:) = interp1(standardLine(realCols),turbidity(ii,realCols),  standardLine, 'nearest','extrap'); 
          elseif length(realCols)==1
            fillProfChlo2(ii,:) = chlorophyll(ii,realCols);   
            fillProfOcon2(ii,:) = oxygen_con(ii,realCols); 
            fillProfOsat2(ii,:) = oxygen_sat(ii,realCols); 
            fillProfTurb2(ii,:) = turbidity(ii,realCols); 
          end
        end
      end

      % Delete to bathy depths across channel
      gvBathyMask = NaN(size(ptemp, 1), size(ptemp,2));
      for j = 1:size(gvBathyMask, 2)  
        mask = find (depthLevels >  (round(interpBathy (j))));
        if length(mask)>1       
          fillProfTemp2(mask,j) = NaN;
          fillProfSalt2(mask,j) = NaN;
          if optionbgc
            fillProfChlo2(mask,j) = NaN;
            fillProfOcon2(mask,j) = NaN;
            fillProfOsat2(mask,j) = NaN;
            fillProfTurb2(mask,j) = NaN;
          end
        end
      end   

      % For final plots 
      clear final
      final.ptemp       = fillProfTemp2;
      final.salinity    = fillProfSalt2;
      if optionbgc
        final.chlorophyll = fillProfChlo2;
        final.oxygen_con  = fillProfOcon2;
        final.oxygen_sat  = fillProfOsat2;
        final.turbidity   = fillProfTurb2;
      end       
 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% Smooth T/S prior to GV calc (moving average)

      % Moving average filter to smooth 
      if  rem((smoothParam/binResolutionHoriz),2)==0
        filterSize = (smoothParam/binResolutionHoriz)+1;
      else
        filterSize = round(smoothParam/binResolutionHoriz);
      end    
      newSmoothSignalPT = NaN(length(depthLevels),length(standardLine));
      newSmoothSignalSA = NaN(length(depthLevels),length(standardLine));
      if optionbgc
        newSmoothSignalCH = NaN(length(depthLevels),length(standardLine));
        newSmoothSignalOC = NaN(length(depthLevels),length(standardLine));
        newSmoothSignalOS = NaN(length(depthLevels),length(standardLine));
        newSmoothSignalTU = NaN(length(depthLevels),length(standardLine));
      end
      for j=1:length(depthLevels)
        origSignalPT = fillProfTemp2(j,:);
        origSignalSA = fillProfSalt2(j,:);
        if optionbgc
          origSignalCH = fillProfChlo2(j,:);
          origSignalOC = fillProfOcon2(j,:);
          origSignalOS = fillProfOsat2(j,:);
          origSignalTU = fillProfTurb2(j,:);
        end
        goodCols     = find(~isnan(origSignalPT)); 
        if ~isempty(goodCols)
          origSignalPT = origSignalPT(goodCols);
          origSignalSA = origSignalSA(goodCols);
          newSmoothSignalPT(j,goodCols) = smooth(origSignalPT,filterSize);
          newSmoothSignalSA(j,goodCols) = smooth(origSignalSA,filterSize);
          if optionbgc
            origSignalCH = origSignalCH(goodCols);
            origSignalOC = origSignalOC(goodCols);
            origSignalOS = origSignalOS(goodCols);
            origSignalTU = origSignalTU(goodCols);
            newSmoothSignalCH(j,goodCols) = smooth(origSignalCH,filterSize);
            newSmoothSignalOC(j,goodCols) = smooth(origSignalOC,filterSize);
            newSmoothSignalOS(j,goodCols) = smooth(origSignalOS,filterSize);
            newSmoothSignalTU(j,goodCols) = smooth(origSignalTU,filterSize);
          end         
        end  
      end
    
      % To tidy up
      clear binProfTemp binProfSalt binProfChlo binProfOcon binProfOsat binProfTurb
      clear interpProfTemp interpProfSalt interpProfChlo interpProfOcon interpProfOsat interpProfTurb
      clear interpProfTemp2 interpProfSalt2 interpProfChlo2 interpProfOcon2 interpProfOsat2 interpProfTurb2
      clear fillProfTemp fillProfSalt fillProfChlo fillProfOcon fillProfOsat fillProfTurb
      clear fillProfSalt2 fillProfTemp2 fillProfChlo2 fillProfOcon2 fillProfOsat2 fillProfTurb2

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% Preprocessing for calculate geostrophic currents 

      % Thermal wind with zero ref common depth for pairs of profiles 
      clear gvMatrix meanDepthLevels distMean densMatrix latMeanGV lonMeanGV
      [gvMatrix, meanDepthLevels, distMean, densMatrix latMeanGV lonMeanGV] = calculate_geostrophy(newSmoothSignalPT, newSmoothSignalSA, depthLevels, standardLine, binResolutionHoriz, binResolutionVert,...
                                                                                       latPts, lonPts, omega, dbar2Pascal, options, useRefVels, matoutDir, missionName );
      % Fill missing first line with 2nd line
      gvMatrix(1,:) = gvMatrix(2, :);

      % Interpolate gv back to original depths 
      interpGV = NaN(length(depthLevels), length(distMean)); 
      for j=1:length(distMean)
        interpGV(:,j) = interp1(meanDepthLevels(1,:), gvMatrix(:,j), depthLevels(1,:), 'linear', NaN);
      end
      clear gvMatrixFill

      % Fill to channel and delete bathy, & choose the grid to save
      % Provides option to keep mean position or interp back to T & S bin positions 
      % Generally select option as optimises the channel width
    
      if options.gvGridOrig  

        % Fill gv to channel based on T/S bin/interp grid positions
        gvMatrixFill = NaN(size(newSmoothSignalPT, 1), size(newSmoothSignalPT,2));
        for ii=1:size(interpGV,1)
          realCols = find(~isnan(interpGV(ii,:)));
          if length(realCols)>=2
            gvMatrixFill(ii,:) = interp1(distMean(realCols),interpGV(ii,realCols),standardLine, 'nearest','extrap');       
          end
        end
        % Delete to bathy depths across channel
        gvBathyMask = NaN(size(newSmoothSignalPT, 1), size(newSmoothSignalPT,2));
        for j = 1:size(gvBathyMask, 2)  
          mask = find (depthLevels >  (round(interpBathy (j))));
          if length(mask)>1
            gvMatrixFill(mask,j) = nan;
          end
        end
        clear gvMatrix
        gvMatrix    = gvMatrixFill;
        gvLongitude = lonBin;
        gvLatitude  = latBin;
        gvLine      = standardLine;  
        gvBathy     = interpBathy;
        gvTimeTrans = interpProfTime;
        gvSal       = newSmoothSignalSA;
        gvTemp      = newSmoothSignalPT;
        if optionbgc
          gvChlo      = newSmoothSignalCH;
          gvOcon      = newSmoothSignalOC;
          gvOsat      = newSmoothSignalOS;
          gvTurb      = newSmoothSignalTU;
        end    

      else

        % Keep to mean positions from gv calc
        gvMatrixFill = NaN(size(interpGV, 1), size(interpGV,2));
        for ii=1:length(interpGV)
          realCols = find(~isnan(interpGV(ii,:)));
          if length(realCols)>=2
            gvMatrixFill(ii,:) = interp1(distMean(realCols),interpGV(ii,realCols),distMean, 'nearest','extrap');
          end
        end
        % Delete to bathy depths across channel
        gvBathyMask   = NaN(size(interpGV, 1), size(interpGV,2));
        interpBathyGV = interp1(standardLine, interpBathy, distMean, 'linear', NaN);
        for j = 1:size(gvBathyMask, 2)  
          mask = find (depthLevels >  (round(interpBathyGV (j))));
          if length(mask)>1
            gvMatrixFill(mask,j) = nan;
          end
        end
        % Need to find mean postions for PT, S, and time 
        gvTemp = NaN(size(gvMatrix, 1), size(gvMatrix,2));
        gvSal  = NaN(size(gvMatrix, 1), size(gvMatrix,2));
        if optionbgc
          gvChlo = NaN(size(gvMatrix, 1), size(gvMatrix,2));
          gvOcon = NaN(size(gvMatrix, 1), size(gvMatrix,2));
          gvOsat = NaN(size(gvMatrix, 1), size(gvMatrix,2));
          gvTurb = NaN(size(gvMatrix, 1), size(gvMatrix,2));
        end
        for ii=1:length(newSmoothSignalPT)
          realCols = find(~isnan(newSmoothSignalPT(ii,:)));
          if length(realCols)>=2
            gvTemp(ii,:) = interp1(standardLine(realCols),newSmoothSignalPT(ii,realCols), distMean, 'linear', NaN);    
            gvSal(ii,:)  = interp1(standardLine(realCols),newSmoothSignalSA(ii,realCols), distMean, 'linear', NaN); 
            if optionbgc
              gvChlo(ii,:) = interp1(standardLine(realCols),newSmoothSignalCH(ii,realCols), distMean, 'linear', NaN);    
              gvOcon(ii,:) = interp1(standardLine(realCols),newSmoothSignalOC(ii,realCols), distMean, 'linear', NaN); 
              gvOsat(ii,:) = interp1(standardLine(realCols),newSmoothSignalOS(ii,realCols), distMean, 'linear', NaN);    
              gvTurb(ii,:) = interp1(standardLine(realCols),newSmoothSignalTU(ii,realCols), distMean, 'linear', NaN); 
            end  
          end
        end
        clear gvMatrix gvLongitude gvLatitude gvLine gvBathy gvTime
        gvMatrix    = gvMatrixFill;
        gvLongitude = lonMeanGV;
        gvLatitude  = latMeanGV;
        gvLine      = distMean;  
        gvBathy     = interpBathyGV;
        gvTimeTrans = interpProfTime;

      end

      % Create new variables
      clear densMatrixGV gvMaxDepth
      gvDensMatrix  = sw_dens(gvSal, gvTemp, 0);
      gvMaxDepth    = NaN(1, length(gvTimeTrans)); 

      for jj= 1:length(gvLine)
        realRows= find (~isnan(gvTemp(:,jj)));
        if ~isempty(realRows)
          gvMaxDepth(1,jj) = depthLevels(max(realRows));
        end
      end

      clear channelInd
      if options.channelTransports  
        channelInd = find(gvMaxDepth >=shelfDepth); % to define shelf / channel 
      else
        channelInd = (1:1:length(gvLongitude));
      end

      % Main ts & gv figs
      if optionplot
        [done] = plot_glider_ts_individual_sections(gvDensMatrix,gvTemp,gvSal,gvTimeTrans,depthLevels,gvLongitude,gvBathy,missionName,transectNo,imageDir,section,datamode,dataset,logo,showplot);
        [done] = plot_glider_density_individual_sections(gvDensMatrix,gvTimeTrans,depthLevels,gvLongitude,gvBathy,missionName,transectNo,imageDir,section,datamode,dataset,logo,showplot);
        [done] = plot_glider_geostrophy_individual_sections(gvDensMatrix,gvMatrix,gvTimeTrans,depthLevels,gvLongitude,gvBathy,missionName,transectNo,imageDir,section,datamode,dataset,logo,showplot);
        if optionbgc
        [done] = plot_glider_bgc_individual_sections(gvDensMatrix,gvChlo,gvOcon,gvOsat,gvTurb,gvTimeTrans,depthLevels,gvLongitude, gvBathy,missionName,transectNo,imageDir,section,datamode,dataset,logo,showplot);
        end
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% Save geostrophic currents 

      % Save data per transect always W to E
      gvLat      = [gvLat gvLatitude] ;
      gvLon      = [gvLon gvLongitude];    
      gvTime     = [gvTime gvTimeTrans];
      gvPtemp    = [gvPtemp gvTemp];
      gvSalt     = [gvSalt gvSal];
      gvDens     = [gvDens gvDensMatrix];
      if optionbgc
      gvChlor = [gvChlor gvChlo];
      gvOconc = [gvOconc gvOcon];
      gvOsatu = [gvOsatu gvOsat];
      gvTurbi = [gvTurbi gvTurb];
      end  
      gv         = [gv gvMatrix];
      gvDepth    = depthLevels;
      gvDist     = [gvDist gvLine];
      gvTransect = [gvTransect (transectNo*(ones(1, size(gvMatrix,2))))];
      tsPtemp    = [tsPtemp final.ptemp];
      tsSalt     = [tsSalt final.salinity];

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% Calculate transports
 
      % Transports 
      binWidth  = binResolutionHoriz * 1000; % in m
      binDepth  = binResolutionVert ;        % in m
      transport = gvMatrix .* (binDepth * binWidth * 10^-6);

      if options.channelTransports  % remove shelf data  
        shelfCols = find(gvMaxDepth < shelfDepth);
        transport(:,shelfCols) = NaN;  
      end

      % Find all +ve and -ve flows
      clear posInd negInd
      posInd = find (gvMatrix > 0 & ~isnan(transport));
      negInd = find (gvMatrix < 0 & ~isnan(transport));
      posTransportTot = sum(transport(posInd));
      negTransportTot = sum(transport(negInd));
      index = {'totN','totS'}; % total transports (North & South)
      transports.totalIndex = index;
      transports.total(count,1) = posTransportTot;
      transports.total(count,2) = negTransportTot;

      % Transports by water mass 
      % Fixed-range: Index WIW < 13 deg (Lopez-Jurado et al, 2008; Pinot & Ganachaud, 1999; Pinot et al, 2002)    
      % Geometric: Juza et al. (2018)
      water    = {'AWo','AWr','LIW', 'WIW', 'WDW'};
      year=str2num(missionName(end-3:end));
      if strcmp(critWIW, 'geometry')
        clear gvDepth2
        gvDepth2=repmat(gvDepth,[size(gvTemp,2),1]);
        indWIW = criterion_geometry_detection_WIW(gvTemp',gvSal',gvDepth2);
        indWIW=(+indWIW)>0;
        % Minimum of profiles with WIW for consistence
        prof_WIW=sum(indWIW,2);
        nbprof_WIW=length(find(prof_WIW>0));
        if nbprof_WIW<5
        indWIW(:,:)=0;
        end
        WIW = find(indWIW');
        AWo = find(~indWIW' & gvSal < 38.35 & gvSal >= 37.5);
        AWr = find(~indWIW' & gvSal < 37.5);
        LIW = find(~indWIW' & gvDepth2' < 750  & gvSal >= 38.35);
        WDW = find(~indWIW' & gvDepth2' >= 750 & gvSal >= 38.35);
      else
        limT=13;
        if strcmp(critWIW, 'visual') 
          if year >= 2015 & year < 2016 
          limT=13.15;
          elseif year >= 2016
          limT=13.25; 
          end
        end
        AWo = find(gvTemp >= limT & gvSal < 38.35 & gvSal >= 37.5) ;
        AWr = find(gvTemp >= limT & gvSal < 37.5);
        LIW = find(gvTemp >= limT & gvSal >= 38.35);
        WIW = find(gvTemp <  limT & gvSal <= 38.45) ;
        WDW = find(gvTemp <  limT & gvSal > 38.45);
      end

      if optionplot
      plot_glider_ts_diagram_watermass(gvTemp,gvSal,gvTimeTrans,AWo,AWr,LIW,WIW,WDW,critWIW,missionName,transectNo, imageDir, section,datamode,dataset,logo,showplot);
      end

      gvDepth3=gvDepth2';
      for jj=1:length(water)
        clear ind posInd negInd posTransportInd negTransportInd
        ind = eval(water{1,jj}); 
        posInd = find(gvMatrix(ind)>0);   
        negInd = find(gvMatrix(ind)<0);     
        posTransportInd = find(~isnan(transport(ind(posInd)))); 
        negTransportInd = find(~isnan(transport(ind(negInd))));
        posTransportWM  = sum(transport(ind(posInd(posTransportInd))));
        negTransportWM  = sum(transport(ind(negInd(negTransportInd))));   
        % Save transports by water mass
        transports.waterIndex = [water, water];
        transports.water(count,jj) = posTransportWM;
        transports.water(count,jj+length(water)) = negTransportWM;
        statistics.waterIndex = [water];
        statistics.waterp05T(count,jj) = prctile(gvTemp(ind),5);
        statistics.waterp17T(count,jj) = prctile(gvTemp(ind),17);
        statistics.waterp50T(count,jj) = prctile(gvTemp(ind),50);
        statistics.waterp83T(count,jj) = prctile(gvTemp(ind),83);
        statistics.waterp95T(count,jj) = prctile(gvTemp(ind),95);
        statistics.waterp05S(count,jj) = prctile(gvSalt(ind),5);
        statistics.waterp17S(count,jj) = prctile(gvSalt(ind),17);
        statistics.waterp50S(count,jj) = prctile(gvSalt(ind),50);
        statistics.waterp83S(count,jj) = prctile(gvSalt(ind),83);
        statistics.waterp95S(count,jj) = prctile(gvSalt(ind),95);
        statistics.waterp05Z(count,jj) = prctile(gvDepth3(ind),5);
        statistics.waterp17Z(count,jj) = prctile(gvDepth3(ind),17);
        statistics.waterp50Z(count,jj) = prctile(gvDepth3(ind),50);
        statistics.waterp83Z(count,jj) = prctile(gvDepth3(ind),83);
        statistics.waterp95Z(count,jj) = prctile(gvDepth3(ind),95);
      end

      clear AWo AWr WIW LIW WDW

      % Save dates
      transports.dateIndex = {'mean', 'length', 'start', 'end'};
      if options.channelTransports 
        clear realcols
        realcols = find(~isnan(gvTimeTrans(channelInd)));
        if ~isempty(realcols)
          transports.date(count,1) = mean(gvTimeTrans(channelInd(realcols)));
          transports.date(count,2) = abs(gvTimeTrans(channelInd(realcols(end))) - gvTimeTrans(channelInd(realcols(1))));    
          if gvTimeTrans(channelInd(realcols(end))) > gvTimeTrans(channelInd(realcols(1)))
            transports.date(count,3) = gvTimeTrans(channelInd(realcols(end)));
            transports.date(count,4) = gvTimeTrans(channelInd(realcols(1)));
          else
            transports.date(count,3) = gvTimeTrans(channelInd(realcols(1)));
            transports.date(count,4) = gvTimeTrans(channelInd(realcols(end)));  
          end
          transports.what = 'channel';
        else
          transports.date(count,1:4)=NaN;
        end
      else
        clear realcols
        realcols = find(~isnan(gvTimeTrans)); 
        if ~isempty(realcols)
          transports.meanDate(count,1)   = mean(gvTimeTrans(realcols));
          transports.lengthDays(count,2) = abs(gvTimeTrans(realcols(end)) - gvTimeTrans(realcols(1))); 
          if gvTimeTrans(realcols(end)) > gvTimeTrans(realcols(1))
            transports.date(count,3) = gvTimeTrans(realcols(end));
            transports.date(count,4) = gvTimeTrans(realcols(1));
          else
            transports.date(count,3) = gvTimeTrans(realcols(1));
            transports.date(count,4) = gvTimeTrans(realcols(end));  
          end
          transports.what = 'transect';
        end
      end

      transports.dist2line(count,1) = mean(dist2Line);
      transports.transect(count,1)  = transectNo; 
      transports.mission            = missionName;

      % Save transport matrix
      gvTransport = [gvTransport, transport];
      count=count+1;

    end

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Saving data...

  % Save all data files   
  dataGV.time         = gvTime;
  dataGV.longitude    = gvLon;
  dataGV.latitude     = gvLat;
  dataGV.depth        = gvDepth;
  dataGV.ptemp        = gvPtemp;
  dataGV.sal          = gvSalt;
  dataGV.ro           = gvDens;
  dataGV.bathy        = gvBathy;
  if optionbgc
  dataGV.chlorophyll  = gvChlor;
  dataGV.oxygen_con   = gvOconc;
  dataGV.oxygen_sat   = gvOsatu;
  dataGV.turbidity    = gvTurbi;
  end
  dataGV.gv           = gv;
  dataGV.transport    = gvTransport;
  dataGV.dist         = gvDist;
  dataGV.transect     = gvTransect;
  dataGV.processInfo  = [binResolutionVert,binResolutionHoriz, smoothParam];
  dataGV.processIndex = ['binResolutionVert ','binResolutionHoriz ', 'smoothingkm'];
  dataGV.tsPtemp      = tsPtemp;
  dataGV.tsSalt       = tsSalt;

  if optionfile

    if exist('dataGV')
    save(datafilename, 'dataGV','-append')
    end

    if exist('transports')
    save(datafilename,'transports','-append')
    end

    if exist('statistics')
    save(datafilename,'statistics','-append')
    end

  end

  % Clear variables for mission and stop deployment processing logging.
  clear dataGV transports statistics
  disp(['   End time: ' , datestr((now), 'yyyy-mm-dd HH:MM:SS')]);
  
end
