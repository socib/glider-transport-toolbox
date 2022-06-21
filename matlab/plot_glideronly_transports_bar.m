function plot_glideronly_transports(missionList,datamode,sectionName,logo,inputDir,imageDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                         %  
% function plot_glideronly_transports_bar                                                 %
%                             (missionList,config,datamode,sectionName,inputDir,imageDir) %
%                                                                                         %
% Purpose:                                                                                %
% - Plot N & S transports for each transect for glider and interpolated model             %
%    (from files created by compute_transport_glider_model.m)                             % 
%                                                                                         %
% Inputs:                                                                                 %
% - missionList: {'canalesJan2011';'canalesFeb2011';'canalesMar2011'}                     %
% - datamode: 'dt' o 'rt'                                                                 %
% - sectionName: 'IbizaChannel','MallorcaChannel'                                         %
% - logo: 'socib','copernicus','nologo'                                                   %
% - inputDir: input directory                                                             %
% - imageDir: image directory                                                             %
%                                                                                         %
% Author: Melanie Juza - SOCIB                                                            %
%         mjuza@socib.es                                                                  %
%                                                                                         %
% Date of creation: Feb-2018 (last modification: 03-Jul-2018)                             %     
%                                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Manual test
%clear all; close all;
%missionList={'canalesJan2011';'canalesFeb2011';'canalesMar2011';'canalesMay2011';'canalesJun2011';...
%             'canalesMar2012';'canalesMay2012';'canalesJul2012';'canalesAug2012';'canalesNov2012';...
%             'canalesJan2013';'canalesMar2013';'canalesMay2013';'canalesJul2013';'canalesSep2013';'canalesNov2013';...
%             'canalesFeb2014';'canalesApr2014';'canalesJul2014';'canalesOct2014';'canalesNov2014';...
%             'canalesJan2015';'canalesApr2015';'canalesJun2015';'canalesAug2015';'canalesOct2015';'canalesNov2015';...
%             'canalesJan2016';'canalesFeb2016';'canalesApr2016';'canalesJul2016';'canalesSep2016';'canalesOct2016';'canalesDec2016';...
%             'canalesMar2017';'canalesMay2017';'canalesJul2017';'canalesSep2017';'canalesNov2017';...
%             'canalesJan2018';'canalesMar2018';'canalesJun2018';'canalesSep2018';...
%             'canalesFeb2019';'canalesApr2019';'canalesMay2019';'canalesJul2019';'canalesSep2019';'canalesNov2019';...
%             'canalesJan2020';'canalesJun2020';'canalesJul2020';'canalesSep2020';'canalesNov2020';...
%             'canalesJan2021';'canalesMay2021'};
%datamode='dt';
%sectionName='MallorcaChannel';            
%logo='socib';       
%inputDir='/LOCALDATA/GLIDERS/Transports/Transports_smooth24_critWIW_geometry_onlyglider';
%imageDir=[inputDir '/timeseries/'];
%addpath('/home/mjuza/TOOLS.DEV/glidertoolbox/');
%addpath('/home/mjuza/TOOLS.DEV/glidertoolbox/matlab');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data: {1} should be the reference dataset (obs)
dataset='glider';

% Plot options
option.glidertotal = true;  % Plot Net/North/South of total flow  for glider         
option.watermasses = true;  % Plot North/South of watermass flows for glider  

% Max time distance between two sections (in days)
maxdist=92; 

% Load logo
if ~strcmp(logo,'nologo')
ima=get_logo(logo);
extfig='.png';
else
extfig='_nologo.png';
end

% Figure style
font=25;

% Water masses
strwater{1}='AWo';   
strwater{2}='AWr';   
strwater{3}='LIW';   
strwater{4}='WIW';  
strwater{5}='WMDW';  
strwater{6}='total'; 
nbwm=length(strwater)-1;

% Colors
for ind=1:length(strwater)
  col{ind}=colorlist(ind);
end
col2{1}=colorlist(2);
col2{2}=colorlist(1);
col2{3}=colorlist(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Read transports

alldates=[];
nbtransect=0;

% Loop on mission
for im = 1:length(missionList)

  missionName = missionList{im};
  display(['--- ',missionName,' ---'])

  rep=['matfiles_',dataset,'_',datamode]; 
  filename= fullfile(inputDir,rep,missionName,[dataset,'_',missionName,'_transportGV_',datamode,'_',sectionName,'.mat']);

  % Read transport file
  if exist(filename)

      clear transports
      display(['Reading ',dataset])
      load(filename);

      % Read transport data
      if exist('transports') & sum(transports.total(:,1))~=0  
     
        dates=double(transports.date(:,1)); 
        if ~isempty('watermasses')
          watermasses=transports.waterIndex(1:5);
        end

        % Section incrementation
        nbsection=nbtransect+length(transports.transect);
        % Glider sections/dates = reference
        refdates=dates; 
        alldates(nbtransect+1:nbsection)=refdates;
        % Initialization
        transportN_wat_all(1,1:nbwm+1,nbtransect+1:nbsection)=NaN;
        transportS_wat_all(1,1:nbwm+1,nbtransect+1:nbsection)=NaN;

        % Loop on section to check dates (to coincide between dataset)
        for irefsection=1:length(refdates) 
          for isection=1:length(transports.transect);
            if (abs(dates(isection)-refdates(irefsection))<1);
              % Water mass transport 
              for iw=1:nbwm 
                transportN_wat_all(1,iw,nbtransect+irefsection)=double(transports.water(isection,iw));
                transportS_wat_all(1,iw,nbtransect+irefsection)=double(transports.water(isection,iw+nbwm));
              end
              % Total transport
              iw=nbwm+1;
              transportN_wat_all(1,iw,nbtransect+irefsection)=double(transports.total(isection,1));
              transportS_wat_all(1,iw,nbtransect+irefsection)=double(transports.total(isection,2));                           
              break
            end       
          end
        end  % End loop on section

      end  
   
  end

  % For section incrementation
  nbtransect=nbsection;
 
  clear dates

  clear refdates

end % End loop on mission

% Remove NaN (in dates) & bad section in glider (due to no detected salinity error)
% Remove sections where not all datasets available
indt=find(alldates'~=0 & isfinite(squeeze(transportN_wat_all(1,nbwm+1,:))) & isfinite(alldates') & squeeze(transportN_wat_all(1,nbwm+1,:)<3) & squeeze(transportS_wat_all(1,nbwm+1,:)>-3)); 
dates2=alldates(indt); clear alldates
transportN(1,1:nbwm+1,1:length(indt))=transportN_wat_all(:,:,indt); 
transportS(1,1:nbwm+1,1:length(indt))=transportS_wat_all(:,:,indt);
clear transport*all

% Transport= [Northward transport, Southward transport, Net transports];
transport2(1,1,1:nbwm+1,1:length(indt))=transportN(:,:,:);
transport2(2,1,1:nbwm+1,1:length(indt))=transportS(:,:,:);
transport2(3,1,1:nbwm+1,1:length(indt))=transportN(:,:,:)+transportS(:,:,:);
clear transportN transportS

% Introduce Nan point between 2 points where time distance is long (for clearer figure)
indpt=find(diff(dates2)>maxdist);
nbpt=length(dates2)+length(indpt);
ipt0=1;ii=0;
for ipt=indpt
 transport(:,:,:,ipt0+ii:ipt+ii)=transport2(:,:,:,ipt0:ipt);
 transport(:,:,:,ipt+ii+1)=NaN;
 dates(ipt0+ii:ipt+ii)=dates2(ipt0:ipt);
 dates(ipt+ii+1)=NaN;
 ipt0=ipt+1;
 ii=ii+1;
end
transport(:,:,:,ipt0+ii:nbpt)=transport2(:,:,:,ipt0:end);
dates(ipt0+ii:nbpt)=dates2(ipt0:end);

% Check in dates: different dates!
for i=1:length(dates)-1
  if dates(i+1)==dates(i);
    dates (i+1)=dates(i)+1;
  end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Conventional (temporal) statistics 
for itrans=1:3

  %%% Total transport
  % MEAN, STD, RMSE, ME, CRMSD, CORR (between dataset and ref dataset)
  tmpfile=transport(itrans,1,nbwm+1,:);
  tmpref =transport(itrans,1,nbwm+1,:);
  meanT(itrans,1)=mean(tmpfile);
  stdT(itrans,1) =std(tmpfile);

  %%% Water mass transport
  for iw=1:nbwm+1
   % CORR (between water mass and total), % of water mass
   tmpwat=transport(itrans,1,iw,:);
   corT_wm(itrans,1,iw)=corrcoef(squeeze(tmpfile(isfinite(tmpfile))),squeeze(tmpwat(isfinite(tmpwat))));
   if itrans~=3
   prcT_wm(itrans,1,iw)=100*sum(tmpwat(isfinite(tmpwat)))/sum(tmpfile(isfinite(tmpfile)));
   end
   clear tmpwat
  end
  clear tmpfile tmpfile

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Text for figures

%%% Total transport

for itrans=1:3
  clear strmean strstd strrmsd strme strcrms strcorr
  strmean='Mean:';
  strstd ='std:';
  strmean=strcat(strmean,num2str(meanT(itrans,1),'%2.2f'),'/');
  strstd =strcat(strstd, num2str(stdT(itrans,1) ,'%2.2f'),'/');
  strstatT{itrans}=[strmean(1:end-1),' ,',strstd(1:end-1)];
end

%%% Water mass transport
for iw=1:nbwm+1
  strstat_wm{iw} =[strwater{iw},': ',num2str(prcT_wm(1,1,iw),'%2.0f'),'/',num2str(prcT_wm(2,1,iw),'%2.0f'),...
                                   '%, R=',num2str(corT_wm(1,1,iw),'%2.2f'),'/',num2str(corT_wm(2,1,iw),'%2.2f'),' (n/s)'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Figures

if strcmp(sectionName,'IbizaChannel')==1
  mint=-2.5;
  maxt=3;
elseif strcmp(sectionName,'MallorcaChannel')==1
  mint=-1;
  maxt=1.5;
end

%%% Glider: net/north/south of total flow

if option.glidertotal

  figure('visible','on'); initfigall(45,30); clf;
  colormap('jet')
  % Northward transports 
  tmp=squeeze(squeeze(transport(1,1,6,:)));
  h = bar(dates (:), tmp','stack'); hold on
  set(h,'LineStyle','none','BarWidth', 2);
  set(get(h(1),'BaseLine'),'LineWidth',2,'LineStyle',':')
  % Southward transports
  tmp=squeeze(squeeze(transport(2,1,6,:)));
  k = bar(dates (:), tmp','stack'); hold on
  set(k,'LineStyle','none','BarWidth', 2);
  set(get(k(1),'BaseLine'),'LineWidth',2,'LineStyle',':')   
  axis([dates(1) dates(end) mint maxt]);grid on;
  ylabel('Total transport / section (Sv)','fontsize',font)
  title([dataset,' in the ',sectionName,' from ',datestr(dates2(1),'mmm-yyyy'),' to ',datestr(dates2(end),'mmm-yyyy')],'fontsize',font)
  set(gca,'fontsize',font)
  datetick('x','mmmyy')
  axis([dates(1) dates(end) mint maxt]);
  plot([dates(1) dates(end)],[0 0],'k--');
  grid on;
  text(dates(10),maxt-(maxt-mint)/15,['Mean(northward/southward/net) = ',num2str(meanT(1,1),'%2.2f'),'/',num2str(meanT(2,1),'%2.2f'),'/',num2str(meanT(3,1),'%2.2f')],'fontsize',font);
  text(dates(10),maxt-(maxt-mint)/8, ['Std(northward/southward/net) = ',num2str(stdT(1,1),'%2.2f'),'/',num2str(stdT(2,1),'%2.2f'),'/',num2str(stdT(3,1),'%2.2f')],'fontsize',font);
  datetick('keeplimits')
  % Logo
  if ~strcmp(logo,'nologo')
  %axes('position',[0.01,0.9,0.09,0.09])
  %imshow(ima)
  axes('position',[0.035,0.9,0.06,0.09])
  image(ima)  
  h = gca; h.XAxis.Visible = 'off';
  h = gca; h.YAxis.Visible = 'off';
  end
  set(gcf,'renderer','zbuffer');
  % Print
  figurepng=strcat([imageDir,'Transports_',sectionName,'_',datestr(dates2(1),'yyyymm'),'-',datestr(dates2(end),'yyyymm'),'_',dataset,'_',strwater{iw},'_bars',extfig])
  print('-dpng','-r300','-painters',figurepng)

end

%%% Glider north/south of watermass flows
if option.watermasses

  watermasses={'AWo','AWr','LIW','WIW','WMDW'};
  figure('visible','on'); initfigall(45,30); clf;
  colormap('jet')
  % Northward transports 
  tmp=squeeze(squeeze(transport(1,1,1:5,:)));
  h = bar(dates (:), tmp','stack'); hold on
  set(h,'LineStyle','none','BarWidth', 2);
  set(get(h(1),'BaseLine'),'LineWidth',2,'LineStyle',':')
  % Southward transports
  tmp=squeeze(squeeze(transport(2,1,1:5,:)));
  k = bar(dates (:), tmp','stack'); hold on
  set(k,'LineStyle','none','BarWidth', 2);
  set(get(k(1),'BaseLine'),'LineWidth',2,'LineStyle',':')   
  axis([dates(1) dates(end) mint maxt]);grid on;
  ylabel('Transport / section (Sv)','fontsize',font)
  title([dataset,' in the ',sectionName,' from ',datestr(dates2(1),'mmm-yyyy'),' to ',datestr(dates2(end),'mmm-yyyy')],'fontsize',font)
  set(gca,'fontsize',font)
  datetick('x','mmmyy')
  axis([dates(1) dates(end) mint maxt]);
  plot([dates(1) dates(end)],[0 0],'k--');
  legend(watermasses,'Orientation','Horizontal','fontsize',font)
  grid on;
  datetick('keeplimits')
  % Logo
  if ~strcmp(logo,'nologo')
  %axes('position',[0.01,0.9,0.09,0.09])
  %imshow(ima)
  axes('position',[0.035,0.9,0.06,0.09])
  image(ima)  
  h = gca; h.XAxis.Visible = 'off';
  h = gca; h.YAxis.Visible = 'off';
  end
  set(gcf,'renderer','zbuffer');
  % Print
  figurepng=strcat([imageDir,'Transports_',sectionName,'_',datestr(dates2(1),'yyyymm'),'-',datestr(dates2(end),'yyyymm'),'_',dataset,'_watermass_bars',extfig])
  print('-dpng','-r300','-painters',figurepng)

  figure('visible','on'); initfigall(45,30); clf;
  colormap('jet')
  % Northward transports 
  tmp=squeeze(squeeze(transport(1,1,1:5,:)));
  h = bar(dates (:), tmp','stack'); hold on
  set(h,'LineStyle','none','BarWidth', 1.5);
  set(get(h(1),'BaseLine'),'LineWidth',1.5,'LineStyle',':')
  % Southward transports
  tmp=squeeze(squeeze(transport(2,1,1:5,:)));
  k = bar(dates (:), tmp','stack'); hold on
  set(k,'LineStyle','none','BarWidth', 1.5);
  set(get(k(1),'BaseLine'),'LineWidth',1.5,'LineStyle',':')   
  axis([dates(1) dates(end) mint maxt]);grid on;
  ylabel('Transport / section (Sv)','fontsize',font)
  title([dataset,' in the ',sectionName,' from ',datestr(dates2(end)-365,'mmm-yyyy'),' to ',datestr(dates2(end),'mmm-yyyy')],'fontsize',font)
  set(gca,'fontsize',font)
  datetick('x','mmm')
  axis([dates(end)-365 dates(end) mint maxt]);
  plot([dates(end)-365 dates(end)],[0 0],'k--');
  legend(watermasses,'Orientation','Horizontal','fontsize',font)
  grid on;
  datetick('keeplimits')
  % Logo
  if ~strcmp(logo,'nologo')
  %axes('position',[0.01,0.9,0.09,0.09])
  %imshow(ima)
  axes('position',[0.035,0.9,0.06,0.09])
  image(ima)  
  h = gca; h.XAxis.Visible = 'off';
  h = gca; h.YAxis.Visible = 'off';
  end
  set(gcf,'renderer','zbuffer');
  % Print
  figurepng=strcat([imageDir,'Transports_',sectionName,'_',datestr(dates2(1),'yyyymm'),'-',datestr(dates2(end),'yyyymm'),'_',dataset,'_watermass_bars_lastyear',extfig])
  print('-dpng','-r300','-painters',figurepng)

end


   

