function [mode]=PS_DISP_sta(value_type,varargin)

% 28.09.2018	NI	; Calculate mean, standard deviation and precision

% TYPE
% PS_DISP_sta('2d_std') = to calculate and plot mean and standard deviation from 2D mode
% PS_DISP_sta('3d_std') = to calculate and plot mean and standard deviation from 3D mode
% PS_DISp_sta('precision') = to calculate the precision based on the stable area

% Note:
% define "select_location.txt"!

stdargin = nargin ;

% Identify location coordinates
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

if strcmp(value_type,'2d_std')

   load('generate_2d.mat','dU','dE');
   % Find selected scatters using index
   Blat=find(dU(:,2)>min_lat & dU(:,2)<max_lat);
   Blon=find(dU(:,1)>min_lon & dU(:,1)<max_lon);
   Blast=ismember(Blat,Blon);
   indexes=find(Blast);
   index=Blat(indexes);
   X = sprintf('%i scatters have been found on the selected location',length(index));
   disp(X)

   % Create selected data
   for c=1:length(index)
       dU_tes(:,c)=dU(index(c,1),3);
       dE_tes(:,c)=dE(index(c,1),3);
   end
   dU_tes=dU_tes';
   dE_tes=dE_tes';
   % Calculate STD
   mean_dU=mean(dU_tes);
   std_dU=std(dU_tes);
   mean_dE=mean(dE_tes);
   std_dE=std(dE_tes);
   calc_std=[mean_dU;std_dU;mean_dE;std_dE];
   dlmwrite('std.txt',calc_std,'precision',8,'delimiter',' ');
   dlmwrite('index.txt',length(index));

elseif strcmp(value_type,'3d_std')

   load('generate_3d.mat','dU','dE','dN');
   % Find selected scatters using index
   Blat=find(dU(:,2)>min_lat & dU(:,2)<max_lat);
   Blon=find(dU(:,1)>min_lon & dU(:,1)<max_lon);
   Blast=ismember(Blat,Blon);
   indexes=find(Blast);
   index=Blat(indexes);
   X = sprintf('%i scatters have been found on the selected location',length(index));
   disp(X)

   % Create selected data
   for c=1:length(index)
       dU_tes(:,c)=dU(index(c,1),3);
       dE_tes(:,c)=dE(index(c,1),3);
       dN_tes(:,c)=dN(index(c,1),3);
   end
   dU_tes=dU_tes';
   dE_tes=dE_tes';
   dN_tes=dN_tes';
   % Calculate STD
   mean_dU=mean(dU_tes);
   std_dU=std(dU_tes);
   mean_dE=mean(dE_tes);
   std_dE=std(dE_tes);
   mean_dN=mean(dN_tes);
   std_dN=std(dN_tes);
   calc_std=[mean_dU;std_dU;mean_dE;std_dE;mean_dN;std_dN];
   dlmwrite('std.txt',calc_std,'precision',8,'delimiter',' ');
   dlmwrite('index.txt',length(index));

end
