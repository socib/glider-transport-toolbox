function [lonProj, latProj, dist2Prof, gridLine] = stdLineProjection(data, profileIdx, latPts, lonPts, gridResolution)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                %
% function [lonProj, latProj, dist2Prof, gridLine]                               %
%    = stdLineProjection(data, profileIdx, latPts, lonPts, gridResolution)       %
%                                                                                %
% Purpose:                                                                       %
% - This function creates a standard line and project lat/lon onto line          %
%   work out what this is in distance in km along line                           %
%   nominal 0.01 dist used to create line - not too important as will bin later  %
%                                                                                %
% Associated routine:                                                            %
% - projectPointsOnLine.m                                                        %
%                                                                                %
% Authors: SOCIB team (www.socib.es)                                             %                                                           
%                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [yLine,xLine] = interpm(latPts,lonPts,0.01,'lin'); 
 [lonProj, latProj, xFitLine, yFitLine] = projectPointsOnLine(data.longitude(profileIdx)',...
 data.latitude(profileIdx)', xLine, yLine);   

 lineLength = distance([yLine(1,1),yLine(end,1)],[xLine(1,1),xLine(end,1)],'km');
 gridLine = [0:gridResolution:lineLength];
    
 % Convert projected profiles to distances along 'standard' line
 dist2Prof = NaN(1,length(latProj));
 for j=1:length(latProj)  
   dist2Prof(j) = distance([latPts(1,1),latProj(j,1)],[lonPts(1,1),lonProj(j,1)],'km');
 end

end  
    
    
