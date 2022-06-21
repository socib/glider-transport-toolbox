function [interpProfileTemperature, interpProfileSalinity] = interpHoriz(temperature, salinity, binDist, interpMethod, extrapValue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function interpHoriz(temperature, salinity, binDist, interpMethod, extrapValue)
%                                                                     
% Purpose:                                                                           
% - Interpolate variables (temperature, salinity) in horizontal  
%                                                                                      
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-July-2017                  
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpolates using method where 2 profiles, where only one it just leaves
% that profile with no interpolation

interpProfileTemperature = NaN(size(temperature,1), size(binDist,2));
interpProfileSalinity    = interpProfileTemperature;

for i=1:size(temperature,1)

  clear realCols
  realCols = find(~isnan(temperature(i,:)));

  if length(realCols) >1
    interpProfileTemperature(i,:) = interp1(binDist(realCols), temperature(i,realCols), binDist, interpMethod, extrapValue);
    interpProfileSalinity(i,:)    = interp1(binDist(realCols), salinity(i,realCols), binDist, interpMethod, extrapValue);
  elseif length(realCols)==1
    interpProfileTemperature(i,realCols) = temperature(i,realCols);
    interpProfileSalinity(i,realCols)    = salinity(i,realCols);       
  end

end
