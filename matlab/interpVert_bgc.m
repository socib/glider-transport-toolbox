function [interpProfileChlo,interpProfileOcon,interpProfileOsat,interpProfileTurb,interpDepths] = interpVert_bgc(chlorophyll,oxygen_con,oxygen_sat,turbidity,depths,depthResolution,maxDepth,minDepth, interpMethod,extrapValue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function interpHoriz(chlorophyll,oxygen_con,oxygen_sat,turbidity,binDist,interpMethod,extrapValue)
%                                                                     
% Purpose:                                                                           
% - Interpolate variables (chlorophyll,oxygen_con,oxygen_sat,turbidity,depth) in vertical  
%                                                                                      
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-June-2018                  
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depthRange  = (minDepth:depthResolution:maxDepth);
interpProfileChlo = NaN(length(depthRange), size(chlorophyll,2));
interpProfileOcon = interpProfileChlo;
interpProfileOsat = interpProfileChlo;
interpProfileTurb = interpProfileChlo;

for i=1:size(chlorophyll,2)

  clear realRows
  realRows = find(~isnan(oxygen_con(:,i)));
  if length(realRows)>1
  interpProfileChlo(:,i) = interp1(depths(realRows), chlorophyll(realRows,i),depthRange, interpMethod, extrapValue);
  interpProfileOcon(:,i) = interp1(depths(realRows), oxygen_con(realRows,i), depthRange, interpMethod, extrapValue);
  interpProfileOsat(:,i) = interp1(depths(realRows), oxygen_sat(realRows,i), depthRange, interpMethod, extrapValue);
  interpProfileTurb(:,i) = interp1(depths(realRows), turbidity(realRows,i),  depthRange, interpMethod, extrapValue);
  end

end

interpDepths = depthRange;
