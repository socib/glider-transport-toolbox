function [interpProfileChlo,interpProfileOcon,interpProfileOsat,interpProfileTurb] = interpHoriz(chlorophyll,oxygen_con,oxygen_sat,turbidity,binDist,interpMethod,extrapValue)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function interpHoriz(chlorophyll,oxygen_con,oxygen_sat,turbidity, binDist, interpMethod, extrapValue)
%                                                                     
% Purpose:                                                                           
% - Interpolate variables (chlorophyll,oxygen_con,oxygen_sat,turbidity) in horizontal  
%                                                                                                                                                                     
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-June-2018                
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% interpolates using method where 2 profiles, where only one it just leaves
% that profile with no interpolation

interpProfileChlo = NaN(size(chlorophyll,1), size(binDist,2));
interpProfileOcon = interpProfileChlo;
interpProfileOsat = interpProfileChlo;
interpProfileTurb = interpProfileChlo;

for i=1:size(chlorophyll,1)

  clear realCols
  realCols = find(~isnan(oxygen_con(i,:)));

  if length(realCols) >1
    interpProfileChlo(i,:) = interp1(binDist(realCols), chlorophyll(i,realCols),binDist, interpMethod, extrapValue);
    interpProfileOcon(i,:) = interp1(binDist(realCols), oxygen_con(i,realCols), binDist, interpMethod, extrapValue);
    interpProfileOsat(i,:) = interp1(binDist(realCols), oxygen_sat(i,realCols), binDist, interpMethod, extrapValue);
    interpProfileTurb(i,:) = interp1(binDist(realCols), turbidity(i,realCols),  binDist, interpMethod, extrapValue);
  elseif length(realCols)==1
    interpProfileChlo(i,realCols) = chlorophyll(i,realCols);
    interpProfileOcon(i,realCols) = oxygen_con(i,realCols);    
    interpProfileOsat(i,realCols) = oxygen_sat(i,realCols);   
    interpProfileTurb(i,realCols) = turbidity(i,realCols);
  end

end
