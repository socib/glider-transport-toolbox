function [gvMatrix, meanDepthLevels, distMean, densMatrix latMeanGV lonMeanGV] = calculate_geostrophy(ptemp, salinity, depthLevels, standardLine, binResolutionHoriz, depthResolution, latPts, lonPts, omega, dbar2Pascal, options, useRefVels, projectDir, missionName )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function calculate_geostrophy(ptemp, salinity, depthLevels, standardLine, binResolutionHoriz, depthResolution, 
%                 latPts, lonPts, omega, dbar2Pascal, options, useRefVels, projectDir, missionName)
%                                                                     
% Purpose:                                                                           
% - Compute geostrophic velocities for full depths
%   i.e. ref level is a common depth of profile
%   and interpolate all back to orig grid positions                                                                                  
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-July-2017                  
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

densMatrix = sw_dens(salinity, ptemp, 0);

% Find new max depth for data in case of bins and grids
maxDepth = NaN(1, size(ptemp,2)); 
for jj= 1:size(ptemp,2)
  realRows= find (~isnan(ptemp(:,jj))); 
  if ~isempty(realRows)
    maxDepth(1,jj) = depthLevels(max(realRows));
  end
end

% Find mean positions 
distMean = standardLine(1:end-1) + binResolutionHoriz/2;    
az = distance([latPts(1,1);latPts(1,2)],[lonPts(1,1);lonPts(1,2)],'km');
[latMeanGV, lonMeanGV] = reckon(latPts(1,1),lonPts(1,1),km2deg(distMean), az);

% Distance between consecutive profiles (in m) 
L = binResolutionHoriz*1000;   

% f = 2*omega*sin(latitude)   
fcoriolis = 2 * omega * sind(latMeanGV); 

[geopotMatrix, gvMatrix, meanDepthLevels] = calculate_geostrophy_commonDepth(densMatrix, depthLevels, dbar2Pascal, depthResolution, L, fcoriolis, options, useRefVels, projectDir, missionName );
