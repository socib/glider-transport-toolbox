function compute_transports_glider(filename,datamode,datatype,section,smoothParam,critWIW,outputDir,programDir,optionfile,optionplot,optionbgc,logo,showplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    
% function compute_transect_glider(filename,datamode,datatype,section,smoothParam,outputDir,optionplot,showplot) 
%                   
% Inputs:                                                                            
% - filename   = glider file (from SOCIB thredds) [L1 product] 
% - datamode   = 'rt' (real time) o 'dt' (delayed mode)
% - datatype   = 'G' (glider)                                         
% - section    = 'IbizaChannel' or 'MallorcaChannel'
% - smoothParam= smoothing parameter (in km) for moving average smoothing 
% - critWIW    ='fixed-range','visual','geometry'                                                        
% - outputDir  = output directory  
% - optionfile = true or false (to save file)                                     
% - optionplot = true or false (to create figures)
% - logo       = 'socib','copernicus','nologo'
% - showplot   = 'on' or 'off' (to show plot on the screen)
%                                                                     
% Purpose:                                                                           
% - Data processing for transect extraction (call processing_transect_glider.m)                                       
% - Data processing for transport computation (call processing_transport_glider.m) 
%   
% Outputs:
% - matfiles with data interpolated on transect: TS, GV, transport
% - figures: T & S vertical sections, transect map, GV vertical sections
%
%                                                                                 
% Author: Melanie Juza - SOCIB                                                       
%         mjuza@socib.es                                                             
%                                                                                    
% Last modification: 20-Jun-2022)   
%
% History:modifed version v11 (only glider) 
%                                                                                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Manuel test
%clear all;close all; warning off all;
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2011/dep0001_ideep00_ime-sldeep000_L1_2011-01-12_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep01-ime_sldeep001/L1/2011/dep0002_ideep01_ime-sldeep001_L1_2011-02-04_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2011/dep0002_ideep00_ime-sldeep000_L1_2011-02-10_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2011/dep0003_ideep00_ime-sldeep000_L1_2011-03-18_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2011/dep0004_ideep00_ime-sldeep000_L1_2011-05-03_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2011/dep0005_ideep00_ime-sldeep000_L1_2011-06-02_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep02-scb_sgdeep002/L1/2012/dep0001_sdeep02_scb-sgdeep002_L1_2012-03-12_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep02-ime_sldeep002/L1/2012/dep0003_ideep02_ime-sldeep002_L1_2012-05-09_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2012/dep0008_ideep00_ime-sldeep000_L1_2012-07-09_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2012/dep0009_ideep00_ime-sldeep000_L1_2012-08-22_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2012/dep0002_sdeep00_scb-sldeep000_L1_2012-11-27_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2013/dep0003_sdeep00_scb-sldeep000_L1_2013-01-30_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2013/dep0001_sdeep01_scb-sldeep001_L1_2013-03-22_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2013/dep0004_sdeep00_scb-sldeep000_L1_2013-05-20_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2013/dep0005_sdeep00_scb-sldeep000_L1_2013-07-15_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2013/dep0007_sdeep00_scb-sldeep000_L1_2013-09-09_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2013/dep0008_sdeep00_scb-sldeep000_L1_2013-11-01_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2013/dep0009_sdeep00_scb-sldeep000_L1_2013-12-02_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2014/dep0010_sdeep00_scb-sldeep000_L1_2014-02-06_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2014/dep0011_sdeep00_scb-sldeep000_L1_2014-04-07_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2014/dep0012_ideep00_ime-sldeep000_L1_2014-05-25_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2014/dep0012_sdeep00_scb-sldeep000_L1_2014-06-10_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2014/dep0015_sdeep01_scb-sldeep001_L1_2014-07-21_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2014/dep0013_ideep00_ime-sldeep000_L1_2014-10-07_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2014/dep0014_ideep00_ime-sldeep000_L1_2014-11-25_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2015/dep0020_sdeep01_scb-sldeep001_L1_2015-01-28_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2015/dep0014_sdeep00_scb-sldeep000_L1_2015-06-18_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2015/dep0015_sdeep00_scb-sldeep000_L1_2015-08-19_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2015/dep0016_sdeep00_scb-sldeep000_L1_2015-10-19_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/ideep00-ime_sldeep000/L1/2015/dep0015_ideep00_ime-sldeep000_L1_2015-11-03_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2015/dep0017_sdeep00_scb-sldeep000_L1_2015-12-11_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2016/dep0016_sdeep00_scb-sldeep000_L1_2016-01-18_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2016/dep0002_sdeep04_scb-sldeep004_L1_2016-02-23_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2016/dep0018_sdeep00_scb-sldeep000_L1_2016-04-27_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2016/dep0017_sdeep00_scb-sldeep000_L1_2016-05-25_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2016/dep0018_sdeep00_scb-sldeep000_L1_2016-07-12_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2016/dep0003_sdeep04_scb-sldeep004_L1_2016-08-23_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2016/dep0004_sdeep04_scb-sldeep004_L1_2016-09-06_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2016/dep0021_sdeep01_scb-sldeep001_L1_2016-11-04_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2016/dep0022_sdeep01_scb-sldeep001_L1_2016-12-30_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2017/dep0023_sdeep01_scb-sldeep001_L1_2017-02-16_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep05-scb_sldeep005/L1/2017/dep0001_sdeep05_scb-sldeep005_L1_2017-03-03_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2017/dep0024_sdeep01_scb-sldeep001_L1_2017-05-10_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2017/dep0008_sdeep04_scb-sldeep004_L1_2017-07-28_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2017/dep0009_sdeep04_scb-sldeep004_L1_2017-09-13_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep05-scb_sldeep005/L1/2017/dep0002_sdeep05_scb-sldeep005_L1_2017-11-16_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2018/dep0011_sdeep04_scb-sldeep004_L1_2018-01-15_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2018/dep0025_sdeep01_scb-sldeep001_L1_2018-03-19_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2018/dep0021_sdeep00_scb-sldeep000_L1_2018-04-05_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2018/dep0026_sdeep00_scb-sldeep000_L1_2018-06-28_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2019/dep0014_sdeep13_scb-sldeep004_L1_2019-02-04_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep00-scb_sldeep000/L1/2019/dep0027_sdeep00_scb-sldeep000_L1_2019-03-05_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2019/dep0014_sdeep04_scb-sldeep004_L1_2019-04-08_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2019/dep0015_sdeep04_scb-sldeep004_L1_2019-05-15_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2019/dep0016_sdeep04_scb-sldeep004_L1_2019-07-18_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2019/dep0031_sdeep01_scb-sldeep001_L1_2019-09-17_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2019/dep0017_sdeep04_scb-sldeep004_L1_2019-11-07_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2020/dep0032_sdeep01_scb-sldeep001_L1_2020-01-30_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2020/dep0018_sdeep04_scb-sldeep004_L1_2020-06-11_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep01-scb_sldeep001/L1/2020/dep0034_sdeep01_scb-sldeep001_L1_2020-07-29_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2020/dep0019_sdeep04_scb-sldeep004_L1_2020-09-29_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep04-scb_sldeep004/L1/2020/dep0020_sdeep04_scb-sldeep004_L1_2020-11-13_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep06-scb_sldeep006/L1/2021/dep0002_sdeep06_scb-sldeep006_L1_2021-01-12_data_dt.nc';
%filename='http://thredds.socib.es/thredds/dodsC/auv/glider/sdeep06-scb_sldeep006/L1/2021/dep0004_sdeep06_scb-sldeep006_L1_2021-05-17_data_dt.nc';

%datamode='dt';
%datatype='G';
%section='MallorcaChannel';
%smoothParam=24;
%critWIW='geometry';
%outputDir = ['DIRECTORY/Transports_smooth',num2str(smoothParam),'_critWIW_',critWIW,'_onlyglider'];
%programDir=['$HOME/glidertoolbox_JS3_D2PTS'];
%optionbgc = 0;
%optionfile= 1;
%optionplot= 1;
%logo='socib';
%showplot='on';

display(['Computing transport for ',filename])
display('---------------')
display(['Dataset: ',datatype,', ',datamode,', ',section])
display(['Processing parameters: ',num2str(smoothParam),' km smoothing, ',critWIW,' detection for WIW'])
display(['Options: save file (',num2str(optionfile),'), save plots (',num2str(optionplot),'), show plots (',showplot,'), bgc data (',num2str(optionbgc),')'])
display('---------------')

optionfile=logical(optionfile);
optionplot=logical(optionplot);
optionbgc =logical(optionbgc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Path to programs

gliderToolboxDir = fullfile(programDir,'matlab');
addpath(genpath(gliderToolboxDir));
addpath(genpath(programDir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output directory

% Directory 
if datatype == 'G'
  dataset='glider';
end
dirfig=['imagery_',dataset,'_',datamode];
dirmat=['matfiles_',dataset,'_',datamode];

% Make folders for images
matfileDir= fullfile(outputDir,dirmat);
system(['mkdir -p ' matfileDir]);
figureDir = fullfile(outputDir,dirfig);
system(['mkdir -p ' figureDir]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data processing for transport computation

% Quality control processing, vertical interpolation, transect detection 
[datafilename]=processing_transect_glider(filename,datamode,datatype,dataset,section,matfileDir,figureDir,optionbgc,optionfile,optionplot,logo,showplot);
 
% Projection and interpolation into standard line, filling, smoothing
% Computation of geostrophic velocity and transports
[datafilename]=processing_transport_glider(datafilename,datamode,datatype,dataset,section,smoothParam,critWIW,matfileDir,figureDir,optionbgc,optionfile,optionplot,logo,showplot);










