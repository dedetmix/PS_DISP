function [mode]=post_processing(value_type,est_noise,varargin)

% Post Processing Correction
% 20.09.2018	NI	; Correct regional trend and unwrapping error tendency

% TYPE
% post_processing('regional_detrend') = to detrend the regional noise based on stable area
% post_processing('uw_tendency',est_noise) = to correct unwrapping error tendency, default = +/- 5 mm

% Note:
% define "select_location.txt" !

stdargin = nargin ;
load('data_match.mat','var_lonlat','var_vector_asc','var_vector_dsc');
load('interpolate.mat','range');
lambda=0.055465800000000;

if nargin<2
    est_noise=5;
end


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
     % Plot error bar, etc for ascending data
     mean_var_vector_asc=mean(var_vector_asc);
     err=(lambda*1000/4)*ones(size(range));
     diff=zeros(1,length(range));
     for c=2:length(range)
         diff(1,c)=mean_var_vector_asc(1,c)-mean_var_vector_asc(1,c-1);
     end
     diff=diff+mean_var_vector_asc;

     % Correct the absolut difference due to phase unwrapped error | ascending data
     for n=1:length(var_vector_asc)
         for c=2:length(range)
             residu(n,c)=var_vector_asc(n,c)-var_vector_asc(n,c-1);
         end
     end
     for n=1:length(var_vector_asc)
     X = sprintf('Line %i ',n);
     disp(X)
         for c=1:length(range)
             if (residu(n,c) >= (lambda*1000/4)+est_noise) && (residu(1,c) < (lambda*1000/2))
                X = sprintf('number %i date has more than 14 mm displaced, corrected',c);
                disp(X)
                correct(n,c)=residu(n,c)-(lambda*1000/2);
             elseif(residu(n,c) <= (-lambda*1000/4)-est_noise) && (residu(n,c) > (-lambda*1000/2))
                X = sprintf('number %i date has more than -14 mm displaced, corrected',c);
                disp(X)
                correct(n,c)=residu(n,c)+(lambda*1000/2);
             elseif(residu(n,c) >= (lambda*1000/2))
                X = sprintf('number %i date has more than 14x2 mm displaced, indicated BREAK POINT',c);
                disp(X)
                correct(n,c)=residu(n,c)+(lambda*1000);
                disp('changed the value as a break point mark');
                %correct(n,c)=residu(n,c);
             elseif(residu(n,c) <= (-lambda*1000/2))
                X = sprintf('number %i date has more than -14x2 mm displaced, indicated BREAK POINT',c);
                disp(X)
                correct(n,c)=residu(n,c)-(lambda*1000);
                disp('changed the value as a break point mark');
                %correct(n,c)=residu(n,c);
            else
                %X = sprintf('number %i date seems no unwrapping error indication, pass',c);
                %disp(X)
                correct(n,c)=residu(n,c);
            end
         end
     disp(' ')
     end
     % Correct time series displacements
     for n=1:length(var_vector_asc)
         for c=1:length(range)
             if (c==1)
                uw_correct_asc(n,c)=var_vector_asc(n,1)+correct(n,c);
             else
                uw_correct_asc(n,c)=uw_correct_asc(n,c-1)+correct(n,c);
             end
         end
     end

     % Correct the absolut difference due to phase unwrapped error | descending data
     for n=1:length(var_vector_dsc)
         for c=2:length(range)
             residu(n,c)=var_vector_dsc(n,c)-var_vector_dsc(n,c-1);
         end
     end
     for n=1:length(var_vector_dsc)
     X = sprintf('Line %i ',n);
     disp(X)
         for c=1:length(range)
             if (residu(n,c) >= (lambda*1000/4)+est_noise) && (residu(1,c) < (lambda*1000/2))
                X = sprintf('number %i date has more than 14 mm displaced, corrected',c);
                disp(X)
                correct(n,c)=residu(n,c)-(lambda*1000/2);
             elseif(residu(n,c) <= (-lambda*1000/4)-est_noise) && (residu(n,c) > (-lambda*1000/2))
                X = sprintf('number %i date has more than -14 mm displaced, corrected',c);
                disp(X)
                correct(n,c)=residu(n,c)+(lambda*1000/2);
             elseif(residu(n,c) >= (lambda*1000/2))
                X = sprintf('number %i date has more than 14x2 mm displaced, indicated BREAK POINT',c);
                disp(X)
                correct(n,c)=residu(n,c)+(lambda*1000);
                disp('changed the value as a break point mark');
                %correct(n,c)=residu(n,c);
             elseif(residu(n,c) <= (-lambda*1000/2))
                X = sprintf('number %i date has more than -14x2 mm displaced, indicated BREAK POINT',c);
                disp(X)
                correct(n,c)=residu(n,c)-(lambda*1000);
                disp('changed the value as a break point mark');
                %correct(n,c)=residu(n,c);
            else
                %X = sprintf('number %i date seems no unwrapping error indication, pass',c);
                %disp(X)
                correct(n,c)=residu(n,c);
            end
         end
     disp(' ')
     end
     % Correct time series displacements
     for n=1:length(var_vector_dsc)
         for c=1:length(range)
             if (c==1)
                uw_correct_dsc(n,c)=var_vector_dsc(n,1)+correct(n,c);
             else
                uw_correct_dsc(n,c)=uw_correct_dsc(n,c-1)+correct(n,c);
             end
         end
     end
     save('data_match.mat','uw_correct_asc','-append');
     save('data_match.mat','uw_correct_dsc','-append');
end
