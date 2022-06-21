function idx_WIW = criterion_geometry_detection_WIW(temp,salt,dept,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function idx_WIW = criterion_geometry_detection_WIW(temp,salt,dept,varargin)
% 
% This function detects WIW from temperature (T) and salinity (S) profiles 
% using a geometric criterion based on the shape of T/S (Juza et al., 2018)
%
% Inputs :
%   - temp: temperature profiles (first dimension is time, second is depth)
%   - salt: salinity profiles (first dimension is time, second is depth)
%   - dept: depth profiles (first dimension is time, second is depth)
%   - options: 'delta_rho'  -> limits of density
%              'alpha_TS'   -> threshold of the deviation
%              'delta_limS' -> upper limit of salinity
%              'lim_depth'  -> upper limit of depth
%              
% Outputs :
%   - idx_WIW: matrix with 1 where WIW and 0 elsewhere
%
% Authors: Romain Escudier and Melanie Juza (last modifications: June 2018)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
%delta_rho = 0.32;
%alpha_TS = 0.135;
%delta_limS = 0.03;
%lim_depth = 500;

delta_rho = 0.36;
alpha_TS = 0.12;
delta_limS = 0.04;
lim_depth = 500;

% User defined parameters
if ~isempty(varargin)
    id_arg = 1;
    while id_arg <= length(varargin)
        cur_arg = varargin{id_arg};      
        if (strcmpi(cur_arg,'delta_rho')) 
            id_arg = id_arg+1;delta_rho = varargin{id_arg};
        elseif (strcmpi(cur_arg,'alpha_TS')) 
            id_arg = id_arg+1;alpha_TS = varargin{id_arg};          
        elseif (strcmpi(cur_arg,'delta_limS')) 
            id_arg = id_arg+1;delta_limS = varargin{id_arg};
        elseif (strcmpi(cur_arg,'lim_depth')) 
            id_arg = id_arg+1;lim_depth = varargin{id_arg};
        else
            error('Unknown argument : %s\n\n %s',cur_arg,help(mfilename))
        end
        id_arg = id_arg+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrix dimension 
[ntim,ndep] = size(temp);

% Computation of density
dens = sw_dens0(salt,temp);

%%% 1. Determination of LIW density

Tmax=16;
S95 = prctile(salt(:),95);
S99 = prctile(salt(:),99);
tmpD=dens(salt>S95 & salt<S99 & temp<Tmax);
tmpT=temp(salt>S95 & salt<S99 & temp<Tmax);
tmpS=salt(salt>S95 & salt<S99 & temp<Tmax);
dens_LIW = tmpD(find(tmpT==max(tmpT(tmpT<prctile(tmpT(:),99))),1));
clear tmp*

if ~isempty(dens_LIW)

  %%% 2. Identification of deep AW

  % Select area of TS
  lim_dens = [dens_LIW-delta_rho dens_LIW];
  idx_select  = dens>=lim_dens(1) & dens<=lim_dens(2) & dept<lim_depth;
  mask_select = +idx_select; mask_select(~idx_select) = NaN;

  % Compute TS extrema for each profile in the selected area
  I1 = rt_findinmat(idx_select','first');id1 = sub2ind(size(idx_select),(1:ntim),I1);
  I2 = rt_findinmat(idx_select','last') ;id2 = sub2ind(size(idx_select),(1:ntim),I2);
  S1 = salt(id1)';
  S2 = salt(id2)';
  T1 = temp(id1)';
  T2 = temp(id2)';

  % Remove extrema that are too close from each other
  idx_2tooclose = abs(dens(id2)-dens(id1))<0.7*delta_rho;
  S1(idx_2tooclose) = NaN;
  S2(idx_2tooclose) = NaN;
  T1(idx_2tooclose) = NaN;
  T2(idx_2tooclose) = NaN;

  %%% 3. Computation of the line formed by the identified deep AW & LIW

  % Compute line formed by the extrema
  M_Droite = (T2-T1)./(S2-S1);
  M_Droite(M_Droite>3)=NaN;M_Droite(M_Droite<-5)=NaN;
  P_Droite = (T2-M_Droite.*S2);

  %%% 4. Computation of the threshold distance to the line

  % Compute Mindist_TS, function of delta_TS
  delta_TS = sqrt((S2-S1).^2+(T2-T1).^2);
  Mindist_TS = alpha_TS./delta_TS; 

  %%% 5. Determination of WIW

  % Distance to the line > Mindist_TS
  dist_Droite = abs(repmat(M_Droite,[1,ndep]).*salt-temp+repmat(P_Droite,[1,ndep]))./sqrt(1+repmat(M_Droite,[1,ndep]).^2).*mask_select.*(sign(repmat(M_Droite,[1,ndep]).*salt-temp+repmat(P_Droite,[1,ndep])));
  idx_crit2 = dist_Droite > repmat(Mindist_TS,[1,ndep]);

  %%% 6. Extrapolation in formation area

  limT = prctile(temp(idx_crit2),80); 
  limS = [min(salt(idx_crit2))-delta_limS max(salt(idx_crit2))+delta_limS];

  %%% Output

  % WIW indices 
  idx_WIW = +(idx_crit2 & dept<lim_depth);
  if ~isempty(limS) & sum(idx_crit2(:))>10
    idx_WIW((temp<=limT & salt>=limS(1) & salt<=limS(2) & dept<lim_depth)&~idx_WIW & dens>=lim_dens(1) & dens<=lim_dens(2)) = 2; 
  end

else

  idx_WIW=[];

end

