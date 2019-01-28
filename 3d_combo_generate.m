function [ts]=3d_combo_generate(value_type,varargin)

% 28.09.2018	NI	; Generate 3D vectors from LOS ascending, descending and AZI ascending, descending using nearneighbour method
 

% TYPE:
% 3d_combo_generate('disp') to calculate dU,dE,dN at a single different time
% 3d_combo_generate('ts') to calculate dU,dE,dN for time series

stdargin = nargin ;

% NEARNEIGHBOUR method
% prepare the data to have the same size and position (point scatter generalization)
       % import data
       delimiterIn='\t';
       los_asc=importdata('los_asc_nn.xyz',delimiterIn);
       los_dsc=importdata('los_dsc_nn.xyz',delimiterIn);
       azi_asc=importdata('azi_asc_nn.xyz',delimiterIn);
       azi_dsc=importdata('azi_dsc_nn.xyz',delimiterIn);
       if exist('data.mat','file')
          save('data.mat','los_asc','-append');
       else
          save('data.mat','los_asc');
       end
          save('data.mat','los_dsc','-append');
          save('data.mat','azi_asc','-append');
          save('data.mat','azi_dsc','-append');

       % save incedence and azimuth angle
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
       % var_angle: (1)azimuth_asc (2)azimuth_dsc (3)incidence_asc (4)incidence_dsc
       % var_vector: (1)los_asc (2)los_dsc (3)azi_asc (4)azi_dsc
       % var_lonlat: (1)longitude (2)latitude
       i=1;
       for c=1:length(los_asc)
           if (~isnan(los_asc(c,3))) && (~isnan(los_dsc(c,3))) && (~isnan(az_angle_asc(c,3))) && (~isnan(az_angle_dsc(c,3))) && (~isnan(inc_angle_asc(c,3))) && (~isnan(inc_angle_dsc(c,3))) && (~isnan(azi_asc(c,3))) && (~isnan(azi_dsc(c,3))) 
              var_angle(i,:)=[az_angle_asc(c,3),az_angle_dsc(c,3),inc_angle_asc(c,3),inc_angle_dsc(c,3)];
              var_vector(i,:)=[los_asc(c,3),los_dsc(c,3),azi_asc(c,3),azi_dsc(c,3)];
              var_lonlat(i,:)=[az_angle_asc(c,1),az_angle_asc(c,2)];
              i=i+1;
           end
       end
       clear c i;

       if exist('data_match.mat','file')
          save('data_match.mat','var_angle','-append');
       else
          save('data_match.mat','var_angle');
       end
       save('data_match.mat','var_vector','-append');
       save('data_match.mat','var_lonlat','-append');

% Calculate 3D displacement from LOS and AZI data

if strcmp(value_type,'disp')
   X = sprintf('Calculate 3d vectors at a single different time using LOS and AZI');
   disp(X)

   load('data_match.mat')

   for c=1:length(var_vector)
    	% estimate dU,dE,dN with original least square (OLS)
    	A=[var_vector(c,1);var_vector(c,2);var_vector(c,3);var_vector(c,4)];
    	B1=cosd(var_angle(c,3));
	B2=-sind(var_angle(c,3)).*sind(var_angle(c,1)+90); 
	B3=-sind(var_angle(c,3)).*cosd(var_angle(c,1)+90);
	B4=cosd(var_angle(c,4));
	B5=-sind(var_angle(c,4)).*sind(var_angle(c,2)+90);
	B6=-sind(var_angle(c,4)).*cosd(var_angle(c,2)+90);
        B7=0;
        B8=sind(var_angle(c,1)+90);
        B9=cosd(var_angle(c,2)+90);
        B10=0;
        B11=sind(var_angle(c,1)+90);
        B12=cosd(var_angle(c,2)+90);
	B=[B1 B2 B3;B4 B5 B6;B7 B8 B9;B10 B11 B12];
	% calculate m --> [dU;dE;dN] vectors
	[m(:,c),stdx(:,c)]=lscov(B,A);
   end
       
   dU=[var_lonlat(:,1) var_lonlat(:,2) m(1,:)'];
   dE=[var_lonlat(:,1) var_lonlat(:,2) m(2,:)'];
   dN=[var_lonlat(:,1) var_lonlat(:,2) m(3,:)'];
   stdx=stdx';
   clear B1 B2 B3 B4 B5 B6 B7 B8 B9 c A B m;

   % see the vertical scale for plotting
   scale=[min(dU(:,3)) max(dU(:,3))];
   % save data
   dlmwrite('dE.txt',dE,'precision',8,'delimiter',' ');
   dlmwrite('dU.txt',dU,'precision',8,'delimiter',' ');
   dlmwrite('dN.txt',dN,'precision',8,'delimiter',' ');
   dlmwrite('ver_scale.txt',scale,'precision',8,'delimiter',' ');
   if exist('generate_3d_pseudo.mat','file')
	   save('generate_3d_combo.mat','dU','-append');
   else
	   save('generate_3d_combo.mat','dU');
   end
   save('generate_3d_combo.mat','dE','-append');
   save('generate_3d_combo.mat','dN','-append');
   save('generate_3d_combo.mat','stdx','-append');

elseif strcmp(value_type,'ts')
   
   X = sprintf('Calculate 3d vectors for time series using LOS and AZI');
   disp(X)
   disp('Coming Soon ...')

   
end  
