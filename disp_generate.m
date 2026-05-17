% Generate dU,dE,dN with an assumption that dN is parallel to terrain (DEM steepest slope and flow direction)
% 24.01.2018 NI ; start the script
%                 - generate 2D vectors (vertical and eastward)
%                 - generate pseudo 3D vectors (+ northward based on Aspect (slope direction), import from ArcGIS)
% 12.04.2018 NI ; apply for time series data (TS)

% import txt DEM, mag, direction, slope from ArcGIS processing or directly
% do on ArcGIS export format
% delimiterOut = ' ';
% [arah_lonlat,delimiterOut]=importdata('ciloto_mag_dir.txt');
% terrain_dir=[arah_lonlat(:,7),arah_lonlat(:,8),arah_lonlat(:,5)];
% terrain_slope=[arah_lonlat(:,7),arah_lonlat(:,8),arah_lonlat(:,6)];
% dlmwrite('terrain_dir.txt',terrain_dir,'precision',8);
% dlmwrite('terrain_slope.txt',terrain_slope,'precision',8);

% prepare the data to have the same size and position (point scatter generalization)

% (1) LOS ascending and descending
% run on terminal $ prepare_disp3d.sh (mode) 1
delimiterIn='\t';
los_asc=importdata('los_asc_fin.xyz',delimiterIn);
los_dsc=importdata('los_dsc_fin.xyz',delimiterIn);
aspect=importdata('aspect.xyz',delimiterIn);
% dem_az=importdata('dem_az_fin.xyz',delimiterIn);
%delete los_asc_fin.xyz los_dsc_fin.xyz terrain_dir_fin.xyz terrain_slope_fin.xyz
if exist('data.mat','file')
   save('data.mat','los_asc','-append');
else
   save('data.mat','los_asc');
end
save('data.mat','los_dsc','-append');
save('data.mat','aspect','-append');

% (2) Incidence angle both ascending and descending
% run on terminal $ prepare_disp3d.sh (mode) 2
delimiterIn=' ';
inc_angle_asc=importdata('inc_angle_asc.lld',delimiterIn);
az_angle_asc=importdata('az_angle_asc.lld',delimiterIn);
save('data.mat','inc_angle_asc','-append');
save('data.mat','az_angle_asc','-append');
inc_angle_dsc=importdata('inc_angle_dsc.lld',delimiterIn);
az_angle_dsc=importdata('az_angle_dsc.lld',delimiterIn);
save('data.mat','inc_angle_dsc','-append');
save('data.mat','az_angle_dsc','-append');

save('data.mat','mask_re','-append');

% (3) Find the match pixels for each variable
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
% var_angle: (1)azimuth_asc (2)azimuth_dsc (3)incidence_asc (4)incidence_dsc (5)direction
% (6)slope (7)dem_azimuth
% var_vector: (1)los_asc (2)los_dsc
% var_lonlat: (1)longitude (2)latitude
if exist('data_match.mat','file')
   save('data_match.mat','var_angle','-append');
else
   save('data_match.mat','var_angle');
end
save('data_match.mat','var_vector','-append');
save('data_match.mat','var_lonlat','-append');

% (4) Main Computation
for c=1:length(var_vector)
    % step 1--> estimate dU,dE, with original least square (OLS)
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
dlmwrite('dE.txt',dE,'precision',8,'delimiter',' ');
dlmwrite('dU.txt',dU,'precision',8,'delimiter',' ');
if exist('generate.mat','file')
   save('generate.mat','dU','-append');
else
   save('generate.mat','dU');
end
save('generate.mat','dE','-append');

% With look angle formula (Source: Resolving vertical and east-west horizontal motion from differential interferometric synthetic aperture radar The L Aquila earthquake)
% la_asc=90-inc_angle_asc(:,3);
% la_dsc=90-inc_angle_dsc(:,3);
% east=(los_dsc(:,3).*cos(degtorad(la_asc))-(los_asc(:,3).*cos(degtorad(la_dsc))./sin(degtorad(la_asc+la_dsc))));
% vertical=(los_dsc(:,3).*sin(degtorad(la_asc))+(los_asc(:,3).*sin(degtorad(la_dsc))./sin(degtorad(la_asc+la_dsc))));

% plot the 2D vectors using velo_ver.sh and velo_hor.sh on terminal

% (5) Generate north direction based on Aspect Calculation (ArcGIS function, source: DEM)
% define local angle for 4 kuadrant (azimuth to local angle) and calculate dN (psuedo)
clear theta dN dr;
theta=zeros(length(var_angle),1);
dN=zeros(length(var_angle),1);
dr=zeros(length(var_angle),1);
for c=1:length(var_angle)
    if (var_angle(c,5) >= 0) && (var_angle(c,5) <= 90)
	theta(c,1)=90-var_angle(c,5);
	dr(c,1)=abs(dE(c,3)*cosd(theta(c,1)));
	dN(c,1)=dr(c,1)*cosd(90-theta(c,1));
    elseif (var_angle(c,5) > 90) && (var_angle(c,5) <= 180)
	theta(c,1)=var_angle(c,5)-90;
	dr(c,1)=abs(dE(c,3)*cosd(theta(c,1)));
	dN(c,1)=dr(c,1)*cosd(90-theta(c,1));
	dN(c,1)=dN(c,1)*-1;
    elseif (var_angle(c,5) > 180) && (var_angle(c,5) <= 270)
	theta(c,1)=270-var_angle(c,5);
	dr(c,1)=abs(dE(c,3)*cosd(theta(c,1)));
	dN(c,1)=dr(c,1)*cosd(90-theta(c,1));
	dN(c,1)=dN(c,1)*-1;
    elseif (var_angle(c,5) > 270) && (var_angle(c,5) <= 360)
	theta(c,1)=var_angle(c,5)-270;
	dr(c,1)=abs(dE(c,3)*cosd(theta(c,1)));
	dN(c,1)=dr(c,1)*cosd(90-theta(c,1));
    else
	X = sprintf('%i line has a negative azimuth angle, set to quadrant IV',c);
	disp(X)
	theta(c,1)=90+var_angle(c,5);
	dr(c,1)=abs(dE(c,3)*cosd(theta(c,1)));
	dN(c,1)=dr(c,1)*cosd(90-theta(c,1));
    end
    c=c+1;
end
clear c X dr;
dN=[var_lonlat(:,1) var_lonlat(:,2) dN];

% calculate dN (psuedo) (old version)
%dN=dE(:,3)./tan(degtorad(var_angle(:,5)));
%dN=[var_lonlat(:,1) var_lonlat(:,2) dN];
dlmwrite('dN.txt',dN,'precision',8,'delimiter',' ');
save('generate.mat','dN','-append');
save('data_match.mat','theta','-append');
% go to terminal to plot the result $ velo_hor.sh

%%%%%%%%%%%%%%%%%%%%%%%%%%% istirahat 02.02.2018 %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%   mulai 12.04.2018   %%%%%%%%%%%%%%%%%%%%%%%%%

% (1) Interpolate each range of time both asc and dsc to be the same range time
% input ts_v-das_asc.mat & ts_v-das_dsc.mat from STAMPS "ps_plot_ts_v-das.mat"
load('ts_v-das_asc.mat', 'ph_uw', 'day', 'lambda')
ph_disp=-ph_uw(:,:)*lambda*1000/(4*pi);
ts_asc=ph_disp;
day_asc=day;
clear ph_disp ph_uw;
load('ts_v-das_dsc.mat', 'ph_uw', 'day', 'lambda')
ph_disp=-ph_uw(:,:)*lambda*1000/(4*pi);
ts_dsc=ph_disp;
day_dsc=day;
clear ph_disp ph_uw;
range=([min(day_asc):12:max(day_dsc)])'; %depends on which orbit time is beginning, with interval 12 days
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
if exist('interpolate.mat','file')
   save('interpolate.mat','ts_asc','-append');
else
   save('interpolate.mat','ts_asc');
end
save('interpolate.mat','ts_dsc','-append');
save('interpolate.mat','range','-append');
% export to txt
load('ts_v-das_asc.mat', 'lonlat')
dlmwrite('lonlat_asc.txt',lonlat,'precision',8,'delimiter',',');
load('ts_v-das_dsc.mat', 'lonlat')
dlmwrite('lonlat_dsc.txt',lonlat,'precision',8,'delimiter',',');
clear lonlat;
dlmwrite('ts_asc.txt',ts_asc,'precision',8,'delimiter',' ');
dlmwrite('ts_dsc.txt',ts_dsc,'precision',8,'delimiter',' ');

% (2) Compute displacements for time series

% go to terminal, run $ prepare_disp3d.sh 4
% go to matlab, timeseries PATH

delimiterIn=',';
asc_surf=importdata('asc_surface.xyz',delimiterIn);
dsc_surf=importdata('dsc_surface.xyz',delimiterIn);
save('interpolate.mat','asc_surf','-append');
save('interpolate.mat','dsc_surf','-append');

% mask scatters with ADD threshold
load('interpolate.mat', 'asc_surf', 'dsc_surf')
load('../data.mat', 'mask_re', 'az_angle_asc', 'az_angle_dsc', 'inc_angle_asc', 'inc_angle_dsc', 'aspect')
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
clear c i az_angle_asc az_angle_dsc inc_angle_asc inc_angle_dsc aspect;

if exist('data_match.mat','file')
   save('data_match.mat','var_angle','-append');
else
   save('data_match.mat','var_angle');
end
save('data_match.mat','var_vector_asc','-append');
save('data_match.mat','var_vector_dsc','-append');
save('data_match.mat','var_lonlat','-append');

% generate dU,dE & dN
load('interpolate.mat', 'range')
dU_ts=zeros(length(var_angle),length(range));
dE_ts=zeros(length(var_angle),length(range));
dN_ts=zeros(length(var_angle),length(range));

%% dU, dE from InSAR
for n=1:length(range)
    for c=1:length(var_vector_asc)
	% step 1--> estimate dU,dE, with original least square
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

%% dN from aspect
for n=1:length(range)
    for c=1:length(var_angle)
    	if (var_angle(c,5) >= 0) && (var_angle(c,5) <= 90)
	   theta(c,1)=90-var_angle(c,5);
	   dN_ts(c,n)=abs(dE_ts(c,n)*cosd(theta(c,1)))*cosd(90-theta(c,1));
    	elseif (var_angle(c,5) > 90) && (var_angle(c,5) <= 180)
	   theta(c,1)=var_angle(c,5)-90;
	   dN_ts(c,n)=abs(dE_ts(c,n)*cosd(theta(c,1)))*cosd(90-theta(c,1));
	   dN_ts(c,n)=dN_ts(c,n)*-1;
    	elseif (var_angle(c,5) > 180) && (var_angle(c,5) <= 270)
	   theta(c,1)=270-var_angle(c,5);
	   dN_ts(c,n)=abs(dE_ts(c,n)*cosd(theta(c,1)))*cosd(90-theta(c,1));
	   dN_ts(c,n)=dN_ts(c,n)*-1;
    	elseif (var_angle(c,5) > 270) && (var_angle(c,5) <= 360)
	   theta(c,1)=var_angle(c,5)-270;
	   dN_ts(c,n)=abs(dE_ts(c,n)*cosd(theta(c,1)))*cosd(90-theta(c,1));
    	else
	   X = sprintf('time-%i : %i line has a negative azimuth angle, set to quadrant IV',n,c);
	   disp(X)
	   theta(c,1)=90+var_angle(c,5);
	   dN_ts(c,n)=abs(dE_ts(c,n)*cosd(theta(c,1)))*cosd(90-theta(c,1));
    	end
        c=c+1;
    end
    clear c X;
end
clear n;

% adjust the first acq. time to be "0" value
dU_ts_new=zeros(size(dU_ts));
dE_ts_new=zeros(size(dE_ts));
dN_ts_new=zeros(size(dN_ts));
for n=1:length(range)
    dU_ts_new(:,n)=dU_ts(:,n) - dU_ts(:,1);
    dE_ts_new(:,n)=dE_ts(:,n) - dE_ts(:,1);
    dN_ts_new(:,n)=dN_ts(:,n) - dN_ts(:,1);
end
clear n;

date=datetime(range,'ConvertFrom','datenum');
dates=datestr(date);
dlmwrite('date.in',dates);
clear date dates

dlmwrite('dU_ts.txt',dU_ts_new,'precision',8,'delimiter',' ');
dlmwrite('dE_ts.txt',dE_ts_new,'precision',8,'delimiter',' ');
dlmwrite('dN_ts.txt',dN_ts_new,'precision',8,'delimiter',' ');
dlmwrite('lonlat.txt',var_lonlat,'precision',8,'delimiter',' ');
if exist('generate.mat','file')
   save('generate.mat','dU_ts','-append');
else
   save('generate.mat','dU_ts');
end
save('generate.mat','dU_ts_new','-append');
save('generate.mat','dE_ts','-append');
save('generate.mat','dE_ts_new','-append');
save('generate.mat','dN_ts','-append');
save('generate.mat','dN_ts_new','-append');
