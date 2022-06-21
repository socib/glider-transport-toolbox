function [binProfileChlo, binProfileOcon,binProfileOsat, binProfileTurb, binDepths] = binVert_bgc(chlorophyll, oxygen_con, oxygen_sat, turbidity, pressure, depthResolution, maxDepth, minDepth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function binVert(chlorophyll, oxygen_con, oxygen_sat, turbidity, pressure, depthResolution, maxDepth, minDepth)
%                                                                     
% Purpose:                                                                           
% - Bin variables (chlorophyll, oxygen_con, oxygen_sat, turbidity, depth) in vertical  
%                                                                                                                                                                       
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-June-2018                  
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depthLimits = [minDepth, maxDepth];
depthRange  = (minDepth:depthResolution:maxDepth);
depthLevels = length(depthRange);

binDepths   = depthRange+depthResolution/2;
binDepths   = binDepths(1:end-1);
depthLevels = length(binDepths);

topEdge = maxDepth;
botEdge = minDepth;

numBins = floor((topEdge-botEdge)/depthResolution);
binProfileChlo = NaN(numBins, 1);
binProfileOcon = binProfileChlo;
binProfileOsat = binProfileChlo;
binProfileTurb = binProfileChlo;

for i=1:size(chlorophyll,2)

  clear x y1 y2 y3 y4   
  x= pressure;
  y1= chlorophyll(:,i);
  y2= oxygen_con(:,i);
  y3= oxygen_sat(:,i);
  y4= turbidity(:,i);

  binEdges    = linspace(botEdge, topEdge, numBins+1);
  [h,whichBin]= histc(x,binEdges);

  for ii=1:numBins
    flagBinMembers = (whichBin ==ii);
    binMembers1 = y1(flagBinMembers);
    binMembers2 = y2(flagBinMembers);
    binMembers3 = y3(flagBinMembers);
    binMembers4 = y4(flagBinMembers);
    % real
    realRows = ~isnan(y2(flagBinMembers));
    binProfileChlo(ii,i) = mean(binMembers1(realRows));
    binProfileOcon(ii,i) = mean(binMembers2(realRows));
    binProfileOsat(ii,i) = mean(binMembers3(realRows));
    binProfileTurb(ii,i) = mean(binMembers4(realRows));
  end

end

binDepths = depthRange+depthResolution/2;
binDepths = binDepths(1:end-1);


