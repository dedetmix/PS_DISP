function [mode]=post_processing(value_type,varargin)

% Post Processing Correction
% 20.09.2018	NI	; Correct regional trend and unwrapping error tendency

% TYPE
% post_processing('regional_detrend') = to detrend the regional noise based on stable area
% post_processing('uw_tendency') = to correct unwrapping error tendency

% Note:
% define "select_location.txt" !

stdargin = nargin ;
load('data_match.mat','var_lonlat','var_vector_asc','var_vector_dsc');
load('interpolate.mat','range');
lambda=0.055465800000000;

if strcmp(value_type,'select_location')
   % define and collect scatters based on ROI
     fileID = fopen('select_location.txt');
     loc = textscan(fileID,'%s');
     loc = char(loc{1});
     fclose(fileID);

     delimiter='/';
     C = strsplit(loc,delimiter);
     min_lon=str2num(C{1,1});
     max_lon=str2num(C{1,2});
     min_lat=str2num(C{1,3});
     max_lat=str2num(C{1,4});
     clear C loc

     % Find selected scatters using index
     Blat=find(var_lonlat(:,2)>min_lat & var_lonlat(:,2)<max_lat);
     Blon=find(var_lonlat(:,1)>min_lon & var_lonlat(:,1)<max_lon);
     Blast=ismember(Blat,Blon);
     indexes=find(Blast);
     index=Blat(indexes);
     X = sprintf('%i scatters have been found on the selected location',length(index));
     disp(X)

elseif strcmp(value_type,'regional_detrend')
% 1) Removing Noise and Regional Trends (detrend)
     % define location from stabil are with high coherence
     lat_min=-6.728465;
     lat_max=-6.726994;
     lon_min=107.014030;
     lon_max=107.015197;

     % find index
     Blat=find(lonlat(:,2)>lat_min& lonlat(:,2)<lat_max);
     Blon=find(lonlat(:,1)>lon_min & lonlat(:,1)<lon_max);
     Blast=ismember(Blat,Blon);
     indexes=find(Blast);
     index=Blat(indexes);

     % find regional trend from scatters in stabil area
     for c=1:length(index)
         reg_trend_phdisp(:,c)=ph_disp(index(c,1),:);
         reg_trend_phuw(:,c)=ph_uw(index(c,1),:);
     end

     % calculate detrend from ph_uw
     for c=1:length(ph_uw)
         ph_uw_detrend(c,:)=ph_uw(c,:)-mean(reg_trend_phuw);
     end

     % save data

elseif strcmp(value_type,'uw_tendency')
% 2) Detect and Correct Possible Phase Unwrapping Errors
     % Plot error bar, etc
     mean_los_asc_tes=mean(los_asc_tes);
     err=(lambda*1000/4)*ones(size(range));
     diff=zeros(1,length(range));
     for c=2:length(range)
         diff(1,c)=mean_los_asc_tes(1,c)-mean_los_asc_tes(1,c-1);
     end
     diff=diff+mean_los_asc_tes;

     % Correct the absolut difference due to phase unwrapped error
     for c=2:length(range)
         residu(1,c)=mean_los_asc_tes(1,c)-mean_los_asc_tes(1,c-1);
     end
     for c=1:length(range)
         if (residu(1,c) >= (lambda*1000/4)+5) && (residu(1,c) < (lambda*1000/2))
            X = sprintf('number %i date has more than 14 mm displaced, corrected',c);
            disp(X)
            correct(1,c)=residu(1,c)-(lambda*1000/2);
         elseif(residu(1,c) <= (-lambda*1000/4)-5) && (residu(1,c) > (-lambda*1000/2))
            X = sprintf('number %i date has more than -14 mm displaced, corrected',c);
            disp(X)
            correct(1,c)=residu(1,c)+(lambda*1000/2);
         elseif(residu(1,c) >= (lambda*1000/2))
            X = sprintf('number %i date has more than 14x2 mm displaced, indicated BREAK POINT',c);
            disp(X)
            correct(1,c)=residu(1,c)+(lambda*1000);
            disp('changed the value as a break point mark');
         elseif(residu(1,c) <= (-lambda*1000/2))
            X = sprintf('number %i date has more than -14x2 mm displaced, indicated BREAK POINT',c);
            disp(X)
            correct(1,c)=residu(1,c)-(lambda*1000);
            disp('changed the value as a break point mark');
         else
            %X = sprintf('number %i date seems no unwrapping error indication, pass',c);
            %disp(X)
            correct(1,c)=residu(1,c);
         end
     end
     % Correct time series displacements
     for c=1:length(range)
         if (c==1)
            mean_correct(1,c)=mean_los_asc_tes(1,1)+correct(1,c);
         else
            mean_correct(1,c)=mean_correct(1,c-1)+correct(1,c);
         end
     end

end
