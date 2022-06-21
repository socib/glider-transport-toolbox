function [color]=colorlist(ind)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                   %  
% function [color]=colorlist(ind)                                                   %
%                                                                                   %
% Purpose:                                                                          %
% - Define color in this list as function of the index                              %
%                                                                                   %
% Authors: Melanie Juza - Baptiste Mourre - SOCIB                                   %
%         mjuza@socib.es                                                            %
%                                                                                   %
% Date of creation: Jul-2016 (last modification: 06-Jul-2017)                       %     
%                                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steel    = [0.27 0.51 0.71];
burgundy = [0.65 0.16 0.16];
olive    = [0.42 0.56 0.14];
orange   = [1    0.5  0   ];
black    = [0    0    0   ];
grey     = [0.4  0.4  0.4 ];
skyblue  = [0.50 0.55 0.95];
mypink   = [0.9  0.5  0.9 ];
brown    = [0.8  0.4  0.1 ];
purple   = [0.42 0.35 0.80];
red      = [1    0    0   ];

col{1}  = steel;
col{2}  = burgundy;
col{3}  = olive;
col{4}  = orange;
col{5}  = black;
col{6}  = grey;
col{7}  = skyblue;
col{8}  = mypink;
col{9}  = brown;
col{10} = purple;
col{11} = red;

color=col{ind};
