function [binProfileChlo,binProfileOcon,binProfileOsat,binProfileTurb,binProfileTime,binDist] = binHoriz(chlorophyll,oxygen_con,oxygen_sat,turbidity,time,dist2Prof,depthLevels,gridLine,horizResolution)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function binHoriz(chlorophyll,oxygen_con,oxygen_sat,turbudity,time,dist2Prof,depthLevels,gridLine,horizResolution)
%                                                                     
% Purpose:                                                                           
% - Bin variables (chlorophyll,oxygen_con,oxygen_sat,turbudity,time) in horizontal  
%
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-June-2018                  
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Binned data has a mean point resolution
binDist = gridLine + horizResolution/2;
binDist = binDist(1:end-1);

topEdge =   max( gridLine);
botEdge =   min( gridLine);

numBins = (topEdge-botEdge)/horizResolution;
binProfileChlo = NaN(length(depthLevels), numBins);
binProfileOcon = binProfileChlo;
binProfileOsat = binProfileChlo;
binProfileTurb = binProfileChlo;
binProfileTime = NaN(1, numBins);

x = dist2Prof;
y1 = chlorophyll;
y2 = oxygen_con;
y3 = oxygen_sat;
y4 = turbidity;

binEdges = linspace(botEdge, topEdge, numBins+1);
[h,whichBin]= histc(x,binEdges);

for ii=1:numBins
  flagBinMembers = find(whichBin ==ii);
  for jj=1:size(y2,1)
    indReal = find(~isnan(y2(jj,flagBinMembers)));
    if ~isempty(indReal)
      binMembers1 = y1(jj,flagBinMembers(indReal));
      binMembers2 = y2(jj,flagBinMembers(indReal));    
      binMembers3 = y3(jj,flagBinMembers(indReal));
      binMembers4 = y4(jj,flagBinMembers(indReal));  
      binProfileChlo(jj,ii) = mean(binMembers1,2); 
      binProfileOcon(jj,ii) = mean(binMembers2,2);
      binProfileOsat(jj,ii) = mean(binMembers3,2); 
      binProfileTurb(jj,ii) = mean(binMembers4,2);
    end
  end
  binMembers5 = time(1,flagBinMembers);
  binProfileTime(1,ii) = mean(binMembers5,2);
end
