function  [binProfileTemperature, binProfileSalinity, binProfileTime, binDist] = binHoriz(ptemp, salinity, time, dist2Prof ,depthLevels, gridLine, horizResolution)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function binHoriz(ptemp, salinity, time, dist2Prof ,depthLevels, gridLine, horizResolution)
%                                                                     
% Purpose:                                                                           
% - Bin variables (temperature, salinity, time) in horizontal  
%
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-July-2017                  
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Binned data has a mean point resolution
binDist = gridLine + horizResolution/2;
binDist = binDist(1:end-1);

topEdge =   max( gridLine);
botEdge =   min( gridLine);

numBins = (topEdge-botEdge)/horizResolution;
binProfileTemperature   = NaN(length(depthLevels), numBins);
binProfileSalinity   = binProfileTemperature;
binProfileTime = NaN(1, numBins);

x= dist2Prof;
y= ptemp;
z= salinity;

binEdges = linspace(botEdge, topEdge, numBins+1);
[h,whichBin]= histc(x,binEdges);

for ii=1:numBins
  flagBinMembers = find(whichBin ==ii);
  for jj=1:size(y,1)
    indReal = find(~isnan(y(jj,flagBinMembers)));
    if ~isempty(indReal)
      binMembers = y(jj,flagBinMembers(indReal));
      binMembers2 = z(jj,flagBinMembers(indReal));    
      binProfileTemperature(jj,ii) = mean(binMembers,2); 
      binProfileSalinity(jj,ii) = mean(binMembers2,2);
    end
  end
  binMembers3 = time(1,flagBinMembers);
  binProfileTime(1,ii) = mean(binMembers3,2);
end
