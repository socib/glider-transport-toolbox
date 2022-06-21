function [gliderType] = findGliderType(gliderName)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function findGliderType(gliderName)
%                                                                     
% Purpose:                                                                           
% - determine glider type
%                   
% Input:                                                                            
% -  glider name                                                                    
%   
%                                                                               
% Authors: SOCIB team (www.socib.es)                                                                                                                        
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp( gliderName, 'ideep02') || strcmp( gliderName, 'ideep00') || strcmp( gliderName, 'icoast00')
    gliderType = 'g1';
    
elseif strcmp( gliderName, 'sdeep00') || strcmp( gliderName, 'sdeep01') || strcmp( gliderName, 'model')  || strcmp( gliderName, 'sdeep04') || strcmp( gliderName, 'sdeep05')
    gliderType = 'g2';

elseif strcmp( gliderName, 'sdeep02') || strcmp( gliderName, 'sdeep03') 
    gliderType = 'sg'; 

elseif strcmp( gliderName, 'sdeep06') || strcmp( gliderName, 'sdeep07') || strcmp( gliderName, 'sdeep08') || strcmp( gliderName, 'sdeep09')
    gliderType = 'g3';
  
end
