function [asc,dsc]=vslope(value_type,varargin)

% 11.06.2018	NI	; Calculate R-Index and Vslope based on C (the percentage of movement detected along the slope)
%			  Source:
%			  Notti et al, 2014. "A methodology for improving landslide PSI data analysis". International Journal of Remote Sensing
%			  Pizarro et al, 2017. "Mapping Vulnarable urban Areas affected by Slow Moving
%		      	  Landslides". MDPI Remote Sensing

% TYPE:
% vslope('asc') to calculate R-Index, C and Vslope for ascending mean velocity
% vslope('dsc') to calculate R-Index, C and Vslope for descending mean velocity
% vslope('asc_ts_surf') to calculate Vslope for ascending time series Surface
% vslope('dsc_ts_surf') to calculate Vslope for descending time series Surface
% vslope('asc_ts_nn') to calculate Vslope for ascending time series Near-Neighbour
% vslope('dsc_ts_nn') to calculate Vslope for descending time series Near-Neighbour

stdargin = nargin ;

if strcmp(value_type,'asc')
       load('data.mat', 'aspect','slope')
       X = sprintf('V slope from ascending');
       disp(X)

       load('data.mat', 'az_angle_asc','inc_angle_asc', 'los_asc')
       azi=az_angle_asc(:,3);
       inc=inc_angle_asc(:,3);
       los=los_asc;
       clear az_angle_asc inc_angle_asc los_asc;
       los_azi=90+azi;
       % calculate R-Index
       So=sind(aspect(:,3)+270-los_azi);
       Ar=(slope(:,3).*So-inc);
       r_index_asc=-sind(Ar);;

	if exist('V_slope.mat','file')
   		save('V_slope.mat','r_index_asc','-append');
	else
   		save('V_slope.mat','r_index_asc');
	end

	H=cosd(inc);
	N=cosd(90-inc).*cosd(180-los_azi);
	E=cosd(90-inc).*cosd(270-los_azi);

   	C_asc=(cosd(slope(:,3)).*sind(aspect(:,3)-90).*N)+((-1.*cosd(slope(:,3)).*cosd(aspect(:,3)-90)).*E)+(sind(slope(:,3)).*H);
   	vslope_asc=los(:,3)./C_asc;
   	save('V_slope.mat','C_asc','-append');
   	save('V_slope.mat','vslope_asc','-append');

elseif strcmp(value_type,'dsc')
       load('data.mat', 'aspect','slope')
       X = sprintf('V slope from descending');
       disp(X)

       load('data.mat', 'az_angle_dsc','inc_angle_dsc', 'los_dsc')
       azi=az_angle_dsc(:,3);
       inc=inc_angle_dsc(:,3);
       los=los_dsc;
       clear az_angle_dsc inc_angle_dsc los_dsc;
       los_azi=360+(azi+90);
       % calculate R-Index
       So=sind(aspect(:,3)+270-los_azi);
       Ar=(slope(:,3).*So-inc);
       r_index_dsc=-sind(Ar);;

	if exist('V_slope.mat','file')
   		save('V_slope.mat','r_index_dsc','-append');
	else
   		save('V_slope.mat','r_index_dsc');
	end

	H=cosd(inc);
	N=cosd(90-inc).*cosd(180-los_azi);
	E=cosd(90-inc).*cosd(270-los_azi);

   	C_dsc=(cosd(slope(:,3)).*sind(aspect(:,3)-90).*N)+((-1.*cosd(slope(:,3)).*cosd(aspect(:,3)-90)).*E)+(sind(slope(:,3)).*H);
   	vslope_dsc=los(:,3)./C_dsc;
   	save('V_slope.mat','C_dsc','-append');
   	save('V_slope.mat','vslope_dsc','-append')

elseif strcmp(value_type,'asc_ts_surf')
       % go to timeseries PATH, vslope 'asc' must be calculated first
       load('../data.mat', 'aspect','slope')
       X = sprintf('V slope from ascending for time series calculation');
       disp(X)

       load('interpolate.mat', 'asc_surf','range')
       load('../V_slope.mat', 'C_asc')
       asc_ts_surf=zeros(length(C_asc),length(range));
       for c=1:length(range)
       	   asc_ts_surf(:,c)=asc_surf(:,c)./C_asc;
       end
       clear c;

	if exist('V_slope.mat','file')
   		save('V_slope.mat','asc_ts_surf','-append');
	else
   		save('V_slope.mat','asc_ts_surf');
	end

elseif strcmp(value_type,'dsc_ts_surf')
       % go to timeseries PATH, vslope 'dsc' must be calculated first
       load('../data.mat', 'aspect','slope')
       X = sprintf('V slope from descending for time series calculation');
       disp(X)

       load('interpolate.mat', 'dsc_surf','range')
       load('../V_slope.mat', 'C_dsc')
       dsc_ts_surf=zeros(length(C_dsc),length(range));
       for c=1:length(range)
       	   dsc_ts_surf(:,c)=dsc_surf(:,c)./C_dsc;
       end
       clear c;

	if exist('V_slope.mat','file')
   		save('V_slope.mat','dsc_ts_surf','-append');
	else
   		save('V_slope.mat','dsc_ts_surf');
	end

elseif strcmp(value_type,'asc_ts_nn')
       % go to timeseries PATH, vslope 'asc' must be calculated first
       load('../data.mat', 'aspect','slope')
       X = sprintf('V slope from ascending for time series calculation');
       disp(X)

       load('interpolate.mat', 'asc_nn','range')
       load('../V_slope.mat', 'C_asc')
       asc_ts_nn=zeros(length(C_asc),length(range));
       for c=1:length(range)
       	   asc_ts_nn(:,c)=asc_nn(:,c)./C_asc;
       end
       clear c;

	if exist('V_slope.mat','file')
   		save('V_slope.mat','asc_ts_nn','-append');
	else
   		save('V_slope.mat','asc_ts_nn');
	end

elseif strcmp(value_type,'dsc_ts_nn')
       % go to timeseries PATH, vslope 'dsc' must be calculated first
       load('../data.mat', 'aspect','slope')
       X = sprintf('V slope from descending for time series calculation');
       disp(X)

       load('interpolate.mat', 'dsc_nn','range')
       load('../V_slope.mat', 'C_dsc')
       dsc_ts_nn=zeros(length(C_dsc),length(range));
       for c=1:length(range)
       	   dsc_ts_nn(:,c)=dsc_nn(:,c)./C_dsc;
       end
       clear c;

	if exist('V_slope.mat','file')
   		save('V_slope.mat','dsc_ts_nn','-append');
	else
   		save('V_slope.mat','dsc_ts_nn');
	end

else
       X = sprintf('Please specify asc, dsc, asc_ts_surf, dsc_ts_surf, asc_ts_nn or dsc_ts_nn');
       disp(X)
end

% calculate C (percentage)
%scatter(aspect_new,c_asc*100,0.1,'b')
%hold on;
%scatter(aspect_new,c_dsc*100,0.1,'r')
