function [interpProfileTemperature, interpProfileSalinity, interpDepths] = interpVert(temperature, salinity, depths, depthResolution, maxDepth, minDepth, interpMethod, extrapValue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function interpHoriz(temperature, salinity, binDist, interpMethod, extrapValue)
%                                                                     
% Purpose:                                                                           
% - Interpolate variables (temperature, salinity, depth) in vertical  
%                                                                                      
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-July-2017                  
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depthRange  = (minDepth:depthResolution:maxDepth);
interpProfileTemperature = NaN(length(depthRange), size(temperature,2));
interpProfileSalinity    = interpProfileTemperature;

for i=1:size(temperature,2)

  clear realRows
  realRows = find(~isnan(temperature(:,i)));

  interpProfileTemperature(:,i) = interp1(depths(realRows), temperature(realRows,i), depthRange, interpMethod, extrapValue);
  interpProfileSalinity(:,i)    = interp1(depths(realRows), salinity(realRows,i), depthRange, interpMethod, extrapValue);

end

interpDepths = depthRange;
