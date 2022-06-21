function [segmentLon, segmentLat, segmentDepth] = get_bathymetry_sections(latPts, lonPts, dataType, section)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function get_bathymetry_sections(latPts, lonPts, dataType, section)
%                                                                     
% Purpose:                                                                           
% - Read bathymetry and interpolate to standard line   
%                   
% Inputs:                                                                            
% - latitudes 
% - longitudes                                       
% - datatype  = 'G' (glider) or 'M' (model)                                              
% - section = 'IbizaChannel'                                                                     
%   
%                                                                                 
% Authors: SOCIB team (www.socib.es)                                                                                                                         
%                                                                                    
% Last modification: 06-July-2017 (Melanie Juza, mjuza@socib.es)                       
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bathymetry file
%fn ='http://thredds.socib.es/thredds/fileServer/ancillary_data/bathymetry/MED_smith_sandwell_v9_1.nc';
%fn='http://thredds.socib.es/thredds/dodsC/ancillary_data/bathymetry/MED_smith_sandwell_v9_1.nc';
fn='/LOCALDATA/BATHY/DATA/MED_smith_sandwell_v9_1.nc';

%Matlab code example   
lat1 = latPts (1,1);
lat2 = latPts (1,2);
lon1 = lonPts (1,1);
lon2 = lonPts (1,2);

% Assume we have lon1, lon2, lat1 and lat2
lonMin = min(lon1, lon2);
lonMax = max(lon1, lon2);
latMin = min(lat1, lat2);
latMax = max(lat1, lat2);

% Latitude
lat = nc_varget(fn, 'lat');
latRangeIdx = find(lat>= latMin & lat <= latMax);
if isempty(latRangeIdx)
latMin = latMin-0.005;
latMax = latMax+0.005;
latRangeIdx = find(lat>= latMin & lat <= latMax);
end
lat = lat(latRangeIdx);

% Longitude
lon = nc_varget(fn, 'lon');
lonRangeIdx = find(lon >= lonMin & lon <= lonMax);
if isempty(lonRangeIdx) 
lonMin = lonMin-0.005;
lonMax = lonMax+0.005;
lonRangeIdx = find(lon>= lonMin & lat <= lonMax);
end
lon = lon(lonRangeIdx);

% Bathymetry
start = [min(latRangeIdx) , min(lonRangeIdx) ] - ones(1, 2);
count = [length(latRangeIdx), length(lonRangeIdx) ];
bathy = nc_varget(fn, 'topo', start, count);

% Once you have a 2D matrix with bathymetry, you
% can generate the depth profile from there. 
numPoints = max(length(lonRangeIdx), length(latRangeIdx));
[distOverGround, segmentLon, segmentLat] = m_lldist([lon1; lon2], [lat1, lat2], numPoints);
[lonGrid, latGrid] = meshgrid(lon, lat); 
if size(lat,1)>1
    segmentDepth = interp2(lonGrid, latGrid, bathy, segmentLon, segmentLat);
else
    segmentDepth = interp1(lonGrid, double(bathy), segmentLon);
end

end


