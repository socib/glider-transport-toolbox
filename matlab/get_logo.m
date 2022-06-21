function im=get_logo(institution)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_logo load the logo of institution given in aergument
%
% insitution: copernicus 0 socib
%
% Author: Melanie Juza- SOCIB
%         mjuza@socib.es

% Date of creation: 29-Sep-2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp ( institution, 'copernicus') 
  file_logo='/home/mjuza/TOOLS.DEV/glidertoolbox/logo/copernicus-logo.png';
elseif strcmp ( institution, 'socib') 
  file_logo='/home/mjuza/TOOLS.DEV/glidertoolbox/logo/socib-logo.jpg';
end
im=imread(file_logo);

end
