function [mode]=PS_DISP_matlab(value_type,varargin)

% 21.08.2018	NI	; Generate pseudo 2D vectors from los ascending & descending

% TYPE:
% PS_DISP_matlab(number)
% number:
% (1) import gridding data to matlab format
% (2) import azimuth & incedence angle data and select PS scatters match based on Amp. Diff. Dispersion (ADD)
%     for surface method
% (3) import azimuth & incedence data for nearneighbour method
% (4) calculate 2d displacement (vertical & west-eastward) components for the mean velocity with OLS
% (5) interpolate asc & dsc data in time to prepare timeseries calculation
% (6) calculate 2d displacement (vertical & west-eastward) components for timeseries using surface and OLS inverse method
% (7) calculate 2d displacement (vertical & west-eastward) components for timeseries using nearneighbour and OLS inverse method
% (..) 3D displacement vectord are calculated from pseudo_disp_generate.m

stdargin = nargin ;

if strcmp(value_type,'1')

   fileID = fopen('asc.txt');
   asc = textscan(fileID,'%s');
   asc = char(asc{1});
   fclose(fileID);
   load(asc, 'ph_disp')
   fileID2 = fopen('loc_asc.txt');
   loc_asc = textscan(fileID2,'%s');
   loc_asc = char(loc_asc{1});
   fclose(fileID2);
   load(loc_asc, 'lonlat')
   los_asc=[lonlat(:,1) lonlat(:,2) ph_disp];
   dlmwrite('los_asc.txt',los_asc,'precision',8);
   clear ph_disp lonlat fileID fileID2 asc loc_asc;

   fileID = fopen('dsc.txt');
   dsc = textscan(fileID,'%s');
   dsc = char(dsc{1});
   fclose(fileID);
   load(dsc, 'ph_disp')
   fileID2 = fopen('loc_dsc.txt');
   loc_dsc = textscan(fileID2,'%s');
   loc_dsc = char(loc_dsc{1});
   fclose(fileID2);
   load(loc_dsc, 'lonlat')
   los_dsc=[lonlat(:,1) lonlat(:,2) ph_disp];
   dlmwrite('los_dsc.txt',los_dsc,'precision',8);
   clear ph_disp lonlat fileID fileID2 dsc loc_dsc;

elseif strcmp(value_type,'2')

       % SURFACE method
       % prepare the data to have the same size and position (point scatter generalization)
       delimiterIn='\t';
       los_asc=importdata('los_asc_surface.xyz',delimiterIn);
       los_dsc=importdata('los_dsc_surface.xyz',delimiterIn);
       mask_re=importdata('mask_re.xyz',delimiterIn);
       aspect=importdata('aspect.xyz',delimiterIn);
       if exist('data.mat','file')
          save('data.mat','los_asc','-append');
       else
          save('data.mat','los_asc');
       end
          save('data.mat','los_dsc','-append');
          save('data.mat','aspect','-append');
    
       % save azimuth and incedence angle
       delimiterIn='\t';
       inc_angle_asc=importdata('inc_angle_asc.lld',delimiterIn);
       az_angle_asc=importdata('az_angle_asc.lld',delimiterIn);
       save('data.mat','inc_angle_asc','-append');
       save('data.mat','az_angle_asc','-append');
       inc_angle_dsc=importdata('inc_angle_dsc.lld',delimiterIn);
       az_angle_dsc=importdata('az_angle_dsc.lld',delimiterIn);
       save('data.mat','inc_angle_dsc','-append');
       save('data.mat','az_angle_dsc','-append');
       save('data.mat','mask_re','-append');

       % find the match pixels for each variable
       % select coordinates with mask_re has values (not NaN)
       i=1;
       for c=1:length(los_asc)
           if (~isnan(mask_re(c,3)))
               var_angle(i,:)=[az_angle_asc(c,3),az_angle_dsc(c,3),inc_angle_asc(c,3),inc_angle_dsc(c,3),aspect(c,3)];
               var_vector(i,:)=[los_asc(c,3),los_dsc(c,3)];
               var_lonlat(i,:)=[az_angle_asc(c,1),az_angle_asc(c,2)];
               i=i+1;
           end
       end
       clear c i;

       % save the new arangged data
       % var_angle: (1)azimuth_asc (2)azimuth_dsc (3)incidence_asc (4)incidence_dsc (5)aspect
       % var_vector: (1)los_asc (2)los_dsc
       % var_lonlat: (1)longitude (2)latitude
       if exist('data_match.mat','file')
          save('data_match.mat','var_angle','-append');
       else
          save('data_match.mat','var_angle');
       end
       save('data_match.mat','var_vector','-append');
       save('data_match.mat','var_lonlat','-append');

elseif strcmp(value_type,'3')

       % NEARNEIGHBOUR method
       % prepare the data to have the same size and position (point scatter generalization)
       delimiterIn='\t';
       los_asc=importdata('los_asc_nn.xyz',delimiterIn);
       los_dsc=importdata('los_dsc_nn.xyz',delimiterIn);
       aspect=importdata('aspect.xyz',delimiterIn);
       if exist('data.mat','file')
          save('data.mat','los_asc','-append');
       else
          save('data.mat','los_asc');
       end
          save('data.mat','los_dsc','-append');
          save('data.mat','aspect','-append');
    
       % save incedence angle
       delimiterIn='\t';
       inc_angle_asc=importdata('inc_angle_asc.lld',delimiterIn);
       az_angle_asc=importdata('az_angle_asc.lld',delimiterIn);
       save('data.mat','inc_angle_asc','-append');
       save('data.mat','az_angle_asc','-append');
       inc_angle_dsc=importdata('inc_angle_dsc.lld',delimiterIn);
       az_angle_dsc=importdata('az_angle_dsc.lld',delimiterIn);
       save('data.mat','inc_angle_dsc','-append');
       save('data.mat','az_angle_dsc','-append');

       % save the new arangged data
       % var_angle: (1)azimuth_asc (2)azimuth_dsc (3)incidence_asc (4)incidence_dsc (5)aspect
       % var_vector: (1)los_asc (2)los_dsc
       % var_lonlat: (1)longitude (2)latitude
       i=1;
       for c=1:length(los_asc)
           %if (~isnan(los_asc(c,3))) && (~isnan(los_dsc(c,3))) && (~isnan(az_angle_asc(c,3))) && (~isnan(az_angle_dsc(c,3))) && (~isnan(inc_angle_asc(c,3))) && (~isnan(inc_angle_dsc(c,3))) && (~isnan(aspect(c,3))) 
              var_angle(i,:)=[az_angle_asc(c,3),az_angle_dsc(c,3),inc_angle_asc(c,3),inc_angle_dsc(c,3),aspect(c,3)]; % if no aspect file, delete " ,aspect(c,3) " in this line
              var_vector(i,:)=[los_asc(c,3),los_dsc(c,3)];
              var_lonlat(i,:)=[az_angle_asc(c,1),az_angle_asc(c,2)];
               i=i+1;
           %end
       end
       clear c i;

       if exist('data_match.mat','file')
          save('data_match.mat','var_angle','-append');
       else
          save('data_match.mat','var_angle');
       end
       save('data_match.mat','var_vector','-append');
       save('data_match.mat','var_lonlat','-append');

elseif strcmp(value_type,'4')

       % compute 2d displacement (vertical & west-eastward) components for the mean velocity with OLS
       load ('data_match.mat')
       for c=1:length(var_vector)
           % estimate dU,dE, with original least square (OLS)
           A=[var_vector(c,1);var_vector(c,2)];
           B1=cos(degtorad(var_angle(c,3)));
           B2=-sin(degtorad(var_angle(c,3))).*sin(degtorad(var_angle(c,1)+90)); 
                                      % has been checked -> same to 360 - var_angle(c,1) - 270 (Hanssen)
           B3=cos(degtorad(var_angle(c,4)));
           B4=-sin(degtorad(var_angle(c,4))).*sin(degtorad(var_angle(c,2)+90));
           B=[B1 B2;B3 B4];
           % arrange m --> [dU;dE](n)
           m(:,c)=lscov(B,A);
       end
       dU=[var_lonlat(:,1) var_lonlat(:,2) m(1,:)'];
       dE=[var_lonlat(:,1) var_lonlat(:,2) m(2,:)'];
       clear B1 B2 B3 B4 c A B m;

       % see the vertical scale for plotting
       scale=[min(dU(:,3)) max(dU(:,3))];
       % save the 2d result
       dlmwrite('dE.txt',dE,'precision',8,'delimiter',' ');
       dlmwrite('dU.txt',dU,'precision',8,'delimiter',' ');
       dlmwrite('ver_scale.txt',scale,'precision',8,'delimiter',' ');
       if exist('generate_2d.mat','file')
          save('generate_2d.mat','dU','-append');
       else
          save('generate_2d.mat','dU');
       end
       save('generate_2d.mat','dE','-append');

elseif strcmp(value_type,'5')

       file_m = fopen('process.txt');
       process = textscan(file_m, '%s');
       process = char(process{1});
       fclose(file_m); clear file_m;
       if strcmp (process, 'STAMPS')

          % interpolate each range of time both asc and dsc to be the same range time
          % input ts_v-das_asc.mat & ts_v-das_dsc.mat from STAMPS "ps_plot_ts_v-das.mat" , depending to your result, could be V-D / V-DO / V-DA
          % asc
          fileID = fopen('asc_inputTS.txt');
          asc = textscan(fileID,'%s');
          asc = char(asc{1});
          fclose(fileID);
          load(asc, 'ph_uw', 'day', 'lambda', 'lonlat')
          ph_disp=-ph_uw(:,:)*lambda*1000/(4*pi);
          ts_asc=ph_disp;
          day_asc=day;
          dlmwrite('lonlat_asc.txt',lonlat,'precision',8,'delimiter',',');
          clear ph_disp ph_uw lonlat fileID;
          % dsc
          fileID = fopen('dsc_inputTS.txt');
          dsc = textscan(fileID,'%s');
          dsc = char(dsc{1});
          fclose(fileID);
          load(dsc, 'ph_uw', 'day', 'lambda', 'lonlat')
          ph_disp=-ph_uw(:,:)*lambda*1000/(4*pi);
          ts_dsc=ph_disp;
          day_dsc=day;
          dlmwrite('lonlat_dsc.txt',lonlat,'precision',8,'delimiter',',');
          clear ph_disp ph_uw lonlat fileID;
          range=([min(day_asc):12:max(day_dsc)])'; % "min: day_asc & max:day_dsc" depends on which acquisition time beginning, could be otherwise. the interval time : 12 days

       elseif strcmp(process, 'EXTERNAL')

          delimiterIn=',';
          % define and import TS asc file
          ts_asc=importdata('01.txt',delimiterIn);
          % define and import TS dsc file
          ts_dsc=importdata('02.txt',delimiterIn);
          delimiterIn=',';
          % define and import lonlat asc file
          lonlat_asc=importdata('03.txt',delimiterIn);
          % define and import lonlat dsc file
          lonlat_dsc=importdata('04.txt',delimiterIn);
          % define and import asc days file
          dates_asc=importdata('05.txt');
          dates_asc=num2str(dates_asc);
          day_asc=datenum(dates_asc,'yyyymmdd');
          clear dates_asc;
          % define and import dsc days file
          dates_dsc=importdata('06.txt');
          dates_dsc=num2str(dates_dsc);
          day_dsc=datenum(dates_dsc,'yyyymmdd');
          clear dates_dsc;
          % save data
          if exist('ps_plot_ts_external.mat','file')
             save('ps_plot_ts_external.mat','lonlat_asc','-append');
          else
             save('ps_plot_ts_external.mat','lonlat_asc');
          end
             save('ps_plot_ts_external.mat','lonlat_dsc','-append');
             save('ps_plot_ts_external.mat','ts_asc','-append');
             save('ps_plot_ts_external.mat','ts_dsc','-append');
             save('ps_plot_ts_external.mat','day_asc','-append');
             save('ps_plot_ts_external.mat','day_dsc','-append');

       end

       range=([min(day_asc):12:max(day_dsc)])';
       % interpolate ascending
       for c=1:length(ts_asc)
           inter_asc(:,c)=interp1(day_asc,ts_asc(c,:),range,'*linear');
       end
       inter_asc=inter_asc';
       clear c;

       % interpolate descending
       for c=1:length(ts_dsc)
           inter_dsc(:,c)=interp1(day_dsc,ts_dsc(c,:),range,'*linear');
       end
       inter_dsc=inter_dsc';
       clear c;

       % merging, excluding NaN
       i=1;
       for c=1:length(range)
           if (~isnan(inter_asc(1,c))) && (~isnan(inter_dsc(1,c)))
              ts_asc_new(:,i)=[inter_asc(:,c)];
              ts_dsc_new(:,i)=[inter_dsc(:,c)];
              range_new(i,:)=range(c,:);
              i=i+1;
           end
       end
       ts_asc=ts_asc_new;
       ts_dsc=ts_dsc_new;
       range=range_new;
       clear c i ts_asc_new ts_dsc_new day_asc day_dsc inter_asc inter_dsc range_new;

       % save interpolation result
       if exist('interpolate.mat','file')
          save('interpolate.mat','ts_asc','-append');
       else
          save('interpolate.mat','ts_asc');
       end
       save('interpolate.mat','ts_dsc','-append');
       save('interpolate.mat','range','-append');
       % export to txt
       dlmwrite('ts_asc.txt',ts_asc,'precision',8,'delimiter',' ');
       dlmwrite('ts_dsc.txt',ts_dsc,'precision',8,'delimiter',' ');
       clear asc dsc;

elseif strcmp(value_type,'6')

       % calculate 2d displacement (vertical & west-eastward) components for timeseries using surface and OLS inverse method
       delimiterIn=',';
       asc_surf=importdata('asc_surface.xyz',delimiterIn);
       dsc_surf=importdata('dsc_surface.xyz',delimiterIn);
       save('interpolate.mat','asc_surf','-append');
       save('interpolate.mat','dsc_surf','-append');

       % mask scatters with ADD threshold
       load('../data.mat', 'mask_re', 'az_angle_asc', 'az_angle_dsc', 'inc_angle_asc', 'inc_angle_dsc','aspect')
       i=1;
       for c=1:length(asc_surf)
           if (~isnan(mask_re(c,3)))
              var_angle(i,:)=[az_angle_asc(c,3),az_angle_dsc(c,3),inc_angle_asc(c,3),inc_angle_dsc(c,3),aspect(c,3)];
              var_vector_asc(i,:)=[asc_surf(c,:)];
	      var_vector_dsc(i,:)=[dsc_surf(c,:)];
              var_lonlat(i,:)=[az_angle_asc(c,1),az_angle_asc(c,2)];
              i=i+1;
           end
       end
       clear c i az_angle_asc az_angle_dsc inc_angle_asc inc_angle_dsc;

       % save data match for TS
       if exist('data_match.mat','file')
          save('data_match.mat','var_angle','-append');
       else
          save('data_match.mat','var_angle');
       end
       save('data_match.mat','var_vector_asc','-append');
       save('data_match.mat','var_vector_dsc','-append');
       save('data_match.mat','var_lonlat','-append');

       %% generate dU and dE
       load('interpolate.mat', 'range')
       dU_ts=zeros(length(var_angle),length(range));
       dE_ts=zeros(length(var_angle),length(range));

       % dU, dE from InSAR
       for n=1:length(range)
           for c=1:length(var_vector_asc)
	       % estimate dU,dE, with original least square
    	       A=[var_vector_asc(c,n);var_vector_dsc(c,n)];
    	       B1=cos(degtorad(var_angle(c,3)));
    	       B2=-sin(degtorad(var_angle(c,3))).*sin(degtorad(var_angle(c,1)+90));
    	       B3=cos(degtorad(var_angle(c,4)));
    	       B4=-sin(degtorad(var_angle(c,4))).*sin(degtorad(var_angle(c,2)+90));
    	       B=[B1 B2;B3 B4];
   	       % arrange m --> [dU;dE](n)
   	       m(:,c)=lscov(B,A);
           end
           dU_ts(:,n)=[m(1,:)'];
           dE_ts(:,n)=[m(2,:)'];
           clear B1 B2 B3 B4 A B c m;
       end
       clear n;

       % adjust the first acq. time to be "0" value
       dU_ts_new=zeros(size(dU_ts));
       dE_ts_new=zeros(size(dE_ts));
       for n=1:length(range)
           dU_ts_new(:,n)=dU_ts(:,n) - dU_ts(:,1);
           dE_ts_new(:,n)=dE_ts(:,n) - dE_ts(:,1);
       end
       clear n;

       % see the vertical scale for plotting
       Umin=min(dU_ts_new);
       Umax=max(dU_ts_new);
       scale=[min(Umin) max(Umax)];

       % save dates after interpolation
       date=datetime(range,'ConvertFrom','datenum');
       dates=datestr(date);
       dlmwrite('date.in',dates);
       clear date dates

       % export to txt file
       dlmwrite('dU_ts.txt',dU_ts_new,'precision',8,'delimiter',' ');
       dlmwrite('dE_ts.txt',dE_ts_new,'precision',8,'delimiter',' ');
       dlmwrite('lonlat.txt',var_lonlat,'precision',8,'delimiter',' ');
       dlmwrite('ver_scale.txt',scale,'precision',8,'delimiter',' ');
       if exist('generate.mat','file')
          save('generate.mat','dU_ts','-append');
       else
          save('generate.mat','dU_ts');
       end
       save('generate.mat','dU_ts_new','-append');
       save('generate.mat','dE_ts','-append');
       save('generate.mat','dE_ts_new','-append');

elseif strcmp(value_type,'7')

       % calculate 2d displacement (vertical & west-eastward) components for timeseries using nearneighbour and OLS inverse method

       %delimiterIn=',';
       %asc_nn=importdata('asc_nn.xyz',delimiterIn);
       %dsc_nn=importdata('dsc_nn.xyz',delimiterIn);
       load('interpolate.mat', 'range')
       n=length(range)-1;
       format=[repmat('%f,', [1 n]) '%f'];
       asc_nn=cell2mat(textscan(fopen('asc_nn.xyz'),format));
       dsc_nn=cell2mat(textscan(fopen('dsc_nn.xyz'),format));
       clear range n format;
       save('interpolate.mat','asc_nn','-append');
       save('interpolate.mat','dsc_nn','-append');

       % arrange data
       load('../data.mat', 'az_angle_asc', 'az_angle_dsc', 'inc_angle_asc', 'inc_angle_dsc','aspect')
       i=1;
       for c=1:length(asc_nn)
           if (~isnan(asc_nn(c,1))) && (~isnan(dsc_nn(c,1))) && (~isnan(az_angle_asc(c,3))) && (~isnan(az_angle_dsc(c,3))) && (~isnan(inc_angle_asc(c,3))) && (~isnan(inc_angle_dsc(c,3))) && (~isnan(aspect(c,3)))
	      var_angle(i,:)=[az_angle_asc(c,3),az_angle_dsc(c,3),inc_angle_asc(c,3),inc_angle_dsc(c,3),aspect(c,3)];
              var_vector_asc(i,:)=[asc_nn(c,:)];
	      var_vector_dsc(i,:)=[dsc_nn(c,:)];
              var_lonlat(i,:)=[az_angle_asc(c,1),az_angle_asc(c,2)];
              i=i+1;
           end
       end
       clear c i az_angle_asc az_angle_dsc inc_angle_asc inc_angle_dsc;

       % save data match for TS
       if exist('data_match.mat','file')
          save('data_match.mat','var_angle','-append');
       else
          save('data_match.mat','var_angle');
       end
       save('data_match.mat','var_vector_asc','-append');
       save('data_match.mat','var_vector_dsc','-append');
       save('data_match.mat','var_lonlat','-append');

       %% generate dU and dE
       load('interpolate.mat', 'range')
       dU_ts=zeros(length(var_angle),length(range));
       dE_ts=zeros(length(var_angle),length(range));

       % dU, dE from InSAR
       for n=1:length(range)
           for c=1:length(var_vector_asc)
	       % estimate dU,dE, with original least square
    	       A=[var_vector_asc(c,n);var_vector_dsc(c,n)];
    	       B1=cos(degtorad(var_angle(c,3)));
    	       B2=-sin(degtorad(var_angle(c,3))).*sin(degtorad(var_angle(c,1)+90));
    	       B3=cos(degtorad(var_angle(c,4)));
    	       B4=-sin(degtorad(var_angle(c,4))).*sin(degtorad(var_angle(c,2)+90));
    	       B=[B1 B2;B3 B4];
   	       % arrange m --> [dU;dE](n)
   	       m(:,c)=lscov(B,A);
           end
           dU_ts(:,n)=[m(1,:)'];
           dE_ts(:,n)=[m(2,:)'];
           clear B1 B2 B3 B4 A B c m;
       end
       clear n;

       % adjust the first acq. time to be "0" value
       dU_ts_new=zeros(size(dU_ts));
       dE_ts_new=zeros(size(dE_ts));
       for n=1:length(range)
           dU_ts_new(:,n)=dU_ts(:,n) - dU_ts(:,1);
           dE_ts_new(:,n)=dE_ts(:,n) - dE_ts(:,1);
       end
       clear n;

       % save dates after interpolation
       date=datetime(range,'ConvertFrom','datenum');
       dates=datestr(date);
       dlmwrite('date.in',dates);
       clear date dates

       % see the vertical scale for plotting
       Umin=min(dU_ts_new);
       Umax=max(dU_ts_new);
       scale=[min(Umin) max(Umax)];

       % export to txt file
       dlmwrite('dU_ts.txt',dU_ts_new,'precision',8,'delimiter',' ');
       dlmwrite('dE_ts.txt',dE_ts_new,'precision',8,'delimiter',' ');
       dlmwrite('lonlat.txt',var_lonlat,'precision',8,'delimiter',' ');
       dlmwrite('ver_scale.txt',scale,'precision',8,'delimiter',' ');
       if exist('generate.mat','file')
          save('generate.mat','dU_ts','-append');
       else
          save('generate.mat','dU_ts');
       end
       save('generate.mat','dU_ts_new','-append');
       save('generate.mat','dE_ts','-append');
       save('generate.mat','dE_ts_new','-append');

elseif strcmp(value_type,'theta')

       % generate north direction based on Aspect Calculation (ArcGIS function, source: DEM)
       % define local angle for 4 kuadrant (azimuth to local angle)
       load ('data_match.mat')
       theta=zeros(length(var_angle),1);
       for c=1:length(var_angle)
           if (var_angle(c,5) >= 0) && (var_angle(c,5) <= 90)
	      theta(c,1)=90-var_angle(c,5);
           elseif (var_angle(c,5) > 90) && (var_angle(c,5) <= 180)
	      %theta(c,1)=var_angle(c,5)-90;
              theta(c,1)=90-var_angle(c,5);
           elseif (var_angle(c,5) > 180) && (var_angle(c,5) <= 270)
	      %theta(c,1)=270-var_angle(c,5);
              theta(c,1)=90-var_angle(c,5);
           elseif (var_angle(c,5) > 270) && (var_angle(c,5) <= 360)
	      %theta(c,1)=var_angle(c,5)-270;
              theta(c,1)=90-var_angle(c,5);
           elseif (var_angle(c,5) > 360)
	      X = sprintf('%i line has more than 360 azimuth angle, set to quadrant I',c);
	      disp(X)
	      theta(c,1)=90-(var_angle(c,5)-360);
           else
	      X = sprintf('%i line has a negative azimuth angle, set to quadrant IV',c);
	      disp(X)
	      %theta(c,1)=90+var_angle(c,5);
              theta(c,1)=90-(var_angle(c,5)+360);
           end
           c=c+1;
       end
       clear c
       save('data_match.mat','theta','-append');

elseif strcmp(value_type,'azimuth_change') 
       delimiterIn='\t';
       % ascending change
       az_angle_asc=importdata('az_angle_asc.lld',delimiterIn);
       edit_asc=az_angle_asc(:,3)-90;
       az_angle_asc_edit=[az_angle_asc(:,1) az_angle_asc(:,2) edit_asc];
       az_angle_asc=az_angle_asc_edit;
       dlmwrite('az_angle_asc.lld',az_angle_asc,'precision',8,'delimiter','\t');
       % descending change
       az_angle_dsc=importdata('az_angle_dsc.lld',delimiterIn);
       edit_dsc=az_angle_dsc(:,3)-90;
       az_angle_dsc_edit=[az_angle_dsc(:,1) az_angle_dsc(:,2) edit_dsc];
       az_angle_dsc=az_angle_dsc_edit;
       dlmwrite('az_angle_dsc.lld',az_angle_dsc,'precision',8,'delimiter','\t');


end 
