function [lonField, latField, bathyField]  = get_bathymetry_Smith_Sandwell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
% Get bathymetry data from Smith and Sandwell (1997)%
%                                                   %
%  Input arguments:                                 %
%     - filename: bathymetry file                   %
%                                                   %
%  Output arguments:                                %
%     - lonField, latField: 1D field                %
%     - bathyField: 2D field (lat,lon)              %
%                                                   %
% Author: Melanie Juza - SOCIB                      %
%         mjuza@socib.es                            %
% Date of creation: Jan-2016                        %
% Last modification: 07-Jul-2017                    %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fname='http://thredds.socib.es/thredds/fileServer/ancillary_data/bathymetry/MED_smith_sandwell_v9_1.nc';
fname='http://thredds.socib.es/thredds/dodsC/ancillary_data/bathymetry/MED_smith_sandwell_v9_1.nc';
  % Longitude, latitude
  lonField=double(nc_varget(fname,'lon'));
  latField=double(nc_varget(fname,'lat'));

  % Bathymetry
  bathyField=double(nc_varget(fname,'topo')); 


