function [binProfileTemperature, binProfileSalinity, binDepths] = binVert(temperature, salinity, pressure, depthResolution, maxDepth, minDepth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function binVert(temperature, salinity, pressure, depthResolution, maxDepth, minDepth)
%                                                                     
% Purpose:                                                                           
% - Bin variables (temperature, salinity, depth) in vertical  
%                                                                                      
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
% Last modification: 06-July-2017                  
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
binProfileTemperature = NaN(numBins, 1);
binProfileSalinity    = binProfileTemperature;

for i=1:size(temperature,2)

  clear x y z    
  x= pressure;
  y= temperature(:,i);
  z= salinity(:,i);

  binEdges    = linspace(botEdge, topEdge, numBins+1);
  [h,whichBin]= histc(x,binEdges);

  for ii=1:numBins
    flagBinMembers = (whichBin ==ii);
    binMembers = y(flagBinMembers);
    binMembers2 = z(flagBinMembers);
    % real
    realRows = ~isnan(y(flagBinMembers));
    binProfileTemperature(ii,i) = mean(binMembers(realRows));
    binProfileSalinity(ii,i) = mean(binMembers2(realRows));
  end

end

binDepths = depthRange+depthResolution/2;
binDepths = binDepths(1:end-1);
