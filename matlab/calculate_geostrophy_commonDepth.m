function [geopotMatrix, gvMatrix, meanDepthLevels] = calculate_geostrophy_commonDepth(densMatrix, depthDynHeight, dbar2Pascal, interpResolutionVert, L, fcoriolis, options, useRefVels, projectDir, missionName )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function calculate_geostrophy_commonDepth(densMatrix, depthDynHeight, dbar2Pascal, interpResolutionVert, L, fcoriolis, options, useRefVels, projectDir, missionName )
%                                                                     
% Purpose:                                                                           
% - compute the geopotential between pairs of glider profiles in an indexed transect   
%   - deepest common level is the reference level  
%   - sum geopot from min common depth   
%   - find gv                                                                    
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-July-2017                  
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pressure Matrix
pressMatrix = depthDynHeight' * ones(1,size(densMatrix,2)) ; 

% Specific Volume Anomaly: svan = 1/sw_dens(s,t,p) - 1/sw_dens(35,0,p)
svanMatrix = (  ones(size(densMatrix)) ./ densMatrix) - ...
         (ones(size(densMatrix)) ./ sw_dens(35*ones(size(densMatrix)),zeros(size(densMatrix)),pressMatrix) );
     
% Mean anomaly values located at meanDepthLevels
meanSvanMatrix  = svanMatrix((1:end-1),:) + abs(diff(svanMatrix(:,:)))/2; 
meanDepthLevels = depthDynHeight(1:end-1) + interpResolutionVert/2 ;     

% Error assoc. pascal conversion 1.5%
diffPress = diff(dbar2Pascal * pressMatrix); 

% Geopotential
geopot = meanSvanMatrix .* diffPress;        

maxDepthGV   = NaN(1, size(densMatrix, 2));
geopotMatrix = NaN(length(meanDepthLevels), size(densMatrix, 2));
gvMatrix     = NaN(length(meanDepthLevels), size(densMatrix, 2)-1);

for jj= 1:size(densMatrix, 2)
  realRows= find (~isnan(geopot(:,jj)));
  if ~isempty(realRows)
    maxDepthGV(1,jj) = meanDepthLevels(max(realRows));
  end
end

% Need to sum the geopotential from the min common depth
for jj=1:size(densMatrix, 2)-1
 
 clear profile1 profile2 geopotProfile1 geopotProfile2 geopotDiff gv
  commonDepth = min([maxDepthGV(jj) maxDepthGV(jj+1)]); % common depth
  commonRows = find(meanDepthLevels <= commonDepth);    % index of common depth
 
  profile1 =  cumsum(flipud(geopot(commonRows,jj)));   % dynamic ht 
  profile2 =  cumsum(flipud(geopot(commonRows,jj+1))); % dynamic ht 

  geopotProfile1 = flipud(profile1);
  geopotProfile2 = flipud(profile2);

  geopotMatrix(commonRows,jj) = geopotProfile1 ;
  geopotMatrix(commonRows,jj+1) = geopotProfile2;

  geopotDiff = geopotProfile2 - geopotProfile1; 
  gv = geopotDiff  ./ (fcoriolis(jj) .* L);
  gvMatrix(commonRows,jj) =  gv;

end  

if options.useRefVel % not add for shelf
    if useRefVels == 2
     dataName = fullfile(projectDir, [missionName, '_refVel2.mat'])
      load(dataName)
      refLevel = refLevel2
    elseif useRefVels == 1
        refLevel = gvMatrix(101,:) * 0.75
    end
      
      channelInd = find(maxDepthGV >= 200);
      realInd = find(~isnan(refLevel));
      maxDepthGV
      
      if useRefVels == 1
      if channelInd(1) < realInd(1)
          refLevel(channelInd(1):realInd(1)) = refLevel(realInd(1));
      end
      
       if realInd(end) < channelInd(end) 
          refLevel(realInd(end): channelInd(end)) = refLevel(realInd(end));
      end
      end
      
      for pp = 1:length(channelInd)
      profile = channelInd(pp);
      depthVel = refLevel(1,profile);
      vertInd = find(~isnan(gvMatrix(:,profile)));
      gvMatrix(vertInd,profile) = gvMatrix(vertInd,profile)+depthVel;
      end
      
end

end
