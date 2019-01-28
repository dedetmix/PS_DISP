function [mode]=PS_DISP_plot(value_type,varargin)

% 13.09.2018	NI	; Plot 2D/3D/LOS results
% 28.09.2018	NI	; Calculate mean, standard deviation and precision

% TYPE
% PS_DISP_plot('3d_nn') = to plot 3d results from nearneighbour method
% PS_DISP_plot('3d_surf') = to plot 3d results from surface method
% PS_DISP_plot('2d_nn') = to plot 2d results from nearneighbour method
% PS_DISP_plot('2d_surf') = to plot 2d results from surface method
% PS_DISP_plot('los') = to plot LOS asc & dsc results from surface/nearneighbour method
% PS_DISP_plot('2d_std_ts') = to calculate and plot mean and standard deviation from 2D mode for time series
% PS_DISP_plot('3d_std_ts') = to calculate and plot mean and standard deviation from 3D mode for time series

% Note:
% define "select_location.txt" !

stdargin = nargin ;
load('data_match.mat','var_lonlat');
load('interpolate.mat','range');

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

% Find selected scatters using index
Blat=find(var_lonlat(:,2)>min_lat & var_lonlat(:,2)<max_lat);
Blon=find(var_lonlat(:,1)>min_lon & var_lonlat(:,1)<max_lon);
Blast=ismember(Blat,Blon);
indexes=find(Blast);
index=Blat(indexes);
X = sprintf('%i scatters have been found on the selected location',length(index));
disp(X)

% Load data timeseries based on input user
if strcmp(value_type,'3d_nn')

   load('generate_3d_nn.mat','dU_ts_new','dE_ts_new','dN_ts_new');
   % Create selected data
   for c=1:length(index)
       dU_tes(:,c)=dU_ts_new(index(c,1),:);
       dE_tes(:,c)=dE_ts_new(index(c,1),:);
       dN_tes(:,c)=dN_ts_new(index(c,1),:);
   end
   dU_tes=dU_tes';
   dE_tes=dE_tes';
   dN_tes=dN_tes';
   % Plot timeseries scatters
   figure('rend','painters','pos',[10 10 1200 800])
   colormat = rand(length(index),3);
   
   % Plot dU
   subplot(3,1,1);
   for c=1:length(index)
       scatter(range,dU_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dU_tes),'k','LineWidth',1);
   end
   ylabel('dU disp (mm)')
   title('3D results from nearneighbour method')
   %h1 = lsline;
   %h1.LineWidth = 3;
   
   % Plot dE
   subplot(3,1,2);
   for c=1:length(index)
       scatter(range,dE_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dE_tes),'k','LineWidth',1);
   end
   ylabel('dE disp (mm)')
   %h2 = lsline;
   %h2.LineWidth = 3;

   % Plot dN
   subplot(3,1,3);
   for c=1:length(index)
       scatter(range,dN_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dN_tes),'k','LineWidth',1);
   end
   ylabel('dN disp (mm)')
   xlabel('time')
   %h3 = lsline;
   %h3.LineWidth = 3;

   print(gcf,'3d_nn.png','-dpng','-r300');  

elseif strcmp(value_type,'3d_surf')

   load('generate_3d_surface.mat','dU_ts_new','dE_ts_new','dN_ts_new');
   % Create selected data
   for c=1:length(index)
       dU_tes(:,c)=dU_ts_new(index(c,1),:);
       dE_tes(:,c)=dE_ts_new(index(c,1),:);
       dN_tes(:,c)=dN_ts_new(index(c,1),:);
   end
   dU_tes=dU_tes';
   dE_tes=dE_tes';
   dN_tes=dN_tes';
   % Plot timeseries scatters
   figure('rend','painters','pos',[10 10 1200 800])
   colormat = rand(length(index),3);
   
   % Plot dU
   subplot(3,1,1);
   for c=1:length(index)
       scatter(range,dU_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dU_tes),'k','LineWidth',1);
   end
   ylabel('dU disp (mm)')
   title('3D results from surface method')
   %h1 = lsline;
   %h1.LineWidth = 3;
   
   % Plot dE
   subplot(3,1,2);
   for c=1:length(index)
       scatter(range,dE_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dE_tes),'k','LineWidth',1);
   end
   ylabel('dE disp (mm)')
   %h2 = lsline;
   %h2.LineWidth = 3;

   % Plot dN
   subplot(3,1,3);
   for c=1:length(index)
       scatter(range,dN_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dN_tes),'k','LineWidth',1);
   end
   ylabel('dN disp (mm)')
   xlabel('time')
   %h3 = lsline;
   %h3.LineWidth = 3;

   print(gcf,'3d_surf.png','-dpng','-r600');
 
elseif strcmp(value_type,'2d_nn')
   
   load('generate_2d_nn.mat','dU_ts_new','dE_ts_new');
   % Create selected data
   for c=1:length(index)
       dU_tes(:,c)=dU_ts_new(index(c,1),:);
       dE_tes(:,c)=dE_ts_new(index(c,1),:);
   end
   dU_tes=dU_tes';
   dE_tes=dE_tes';
   % Plot timeseries scatters
   figure('rend','painters','pos',[10 10 1200 800])
   colormat = rand(length(index),3);
   
   % Plot dU
   subplot(2,1,1);
   for c=1:length(index)
       scatter(range,dU_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dU_tes),'k','LineWidth',1);
   end
   ylabel('dU disp (mm)')
   title('2D results from nearneighbour method')
   %h1 = lsline;
   %h1.LineWidth = 3;
   
   % Plot dE
   subplot(2,1,2);
   for c=1:length(index)
       scatter(range,dE_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dE_tes),'k','LineWidth',1);
   end
   ylabel('dE disp (mm)')
   xlabel('time')
   %h2 = lsline;
   %h2.LineWidth = 3;

   print(gcf,'2d_nn.png','-dpng','-r600');

elseif strcmp(value_type,'2d_surf')

   load('generate_2d_surface.mat','dU_ts_new','dE_ts_new');
   % Create selected data
   for c=1:length(index)
       dU_tes(:,c)=dU_ts_new(index(c,1),:);
       dE_tes(:,c)=dE_ts_new(index(c,1),:);
   end
   dU_tes=dU_tes';
   dE_tes=dE_tes';
   % Plot timeseries scatters
   figure('rend','painters','pos',[10 10 1200 800])
   colormat = rand(length(index),3);
   
   % Plot dU
   subplot(2,1,1);
   for c=1:length(index)
       scatter(range,dU_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dU_tes),'k','LineWidth',1);
   end
   ylabel('dU disp (mm)')
   title('2D results from surface method')
   %h1 = lsline;
   %h1.LineWidth = 3;
   
   % Plot dE
   subplot(2,1,2);
   for c=1:length(index)
       scatter(range,dE_tes(c,:),5,colormat(c,:),'filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   if (index >= 10)
      plot(range,mean(dE_tes),'k','LineWidth',1);
   end
   ylabel('dE disp (mm)')
   xlabel('time')
   %h2 = lsline;
   %h2.LineWidth = 3;

   print(gcf,'2d_surf.png','-dpng','-r600');

elseif strcmp(value_type,'los')

   % Find any uw_correction or using original var_vector_asc&dsc
   if ~isempty(who('-file', 'data_match.mat', 'uw_correct_asc'))
      load('data_match.mat','uw_correct_asc','uw_correct_dsc');
      var_vector_asc=uw_correct_asc;
      var_vector_dsc=uw_correct_dsc;
      clear uw_correct_asc uw_correct_dsc;
   else
      load('data_match.mat','var_vector_asc','var_vector_dsc');
   end

   % Create selected data
   for c=1:length(index)
       asc_tes(:,c)=var_vector_asc(index(c,1),:);
       dsc_tes(:,c)=var_vector_dsc(index(c,1),:);
   end
   asc_tes=asc_tes';
   dsc_tes=dsc_tes';
   % Plot timeseries scatters
   figure('rend','painters','pos',[10 10 1200 800])
   colormat = rand(length(index),3);
   
   % Plot asc & dsc
   subplot(2,1,1);
   for c=1:length(index)
       scatter(range,asc_tes(c,:),5,'r','filled')
       hold on;
       scatter(range,dsc_tes(c,:),5,'b','filled')
       hold on;
   end
   datetick('x','mm-yy','keepticks')
   ylabel('LOS asc & dsc (mm)')
   title('Displacement from LOS')
   
   % Plot mean asc & dsc
   subplot(2,1,2);
   plot(range,mean(asc_tes),'r','LineWidth',1);
   hold on;
   plot(range,mean(dsc_tes),'b','LineWidth',1);
   hold on;
   scatter(range,mean(asc_tes),8,'k','filled');
   errorbar(range,mean(asc_tes),std(asc_tes),'s','LineWidth',0.1,'MarkerSize',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r')
   hold on;
   scatter(range,mean(dsc_tes),8,'k','filled')
   errorbar(range,mean(dsc_tes),std(dsc_tes),'s','LineWidth',0.1,'MarkerSize',0.5,'MarkerEdgeColor','b','MarkerFaceColor','b')
   hold on;
   lsline
   h=legend('ascending','descending');
   set(h, 'Location', 'southwest')
   datetick('x','mm-yy','keepticks')
   ylabel('mean LOS (mm)')
   xlabel('time')

   print(gcf,'LOS.png','-dpng','-r600');

elseif strcmp(value_type,'2d_std_ts')

   load('generate_2d.mat','dU_ts_new','dE_ts_new');
   % Create selected data
   for c=1:length(index)
       dU_tes(:,c)=dU_ts_new(index(c,1),:);
       dE_tes(:,c)=dE_ts_new(index(c,1),:);
   end
   dU_tes=dU_tes';
   dE_tes=dE_tes';

   % Calculate STD
   mean_dU=mean(dU_tes);
   std_dU=std(dU_tes);
   mean_dE=mean(dE_tes);
   std_dE=std(dE_tes);
   calc_std_2d=[mean_dU;std_dU;mean_dE;std_dE];
   if exist('data.mat','file')
          save('calc_std.mat','calc_std_2d','-append');
   else
          save('calc_std.mat','calc_std_2d');
   end
   dlmwrite('2d_std.txt',calc_std_2d,'precision',8,'delimiter',' ');
   dlmwrite('index.txt',length(index));

   % Plot the mean scatter with error bar
   figure('rend','painters','pos',[10 10 800 600])
   
   % Plot dU
   subplot(2,1,1);
   plot(range,mean_dU,'k','LineWidth',1);
   hold on;
   errorbar(range,mean_dU,std_dU,'s','MarkerSize',1,'MarkerEdgeColor','red','MarkerFaceColor','red')
   datetick('x','mm-yy','keepticks')
   ylabel('dU disp (mm)')
   title('Error bar plot of 2D displacement result')
   
   % Plot dE
   subplot(2,1,2);
   plot(range,mean_dE,'k','LineWidth',1);
   hold on;
   errorbar(range,mean_dE,std_dE,'s','MarkerSize',1,'MarkerEdgeColor','red','MarkerFaceColor','red')
   datetick('x','mm-yy','keepticks')
   ylabel('dE disp (mm)')
   xlabel('time')

   print(gcf,'2d_errorbar.png','-dpng','-r600');

elseif strcmp(value_type,'3d_std_ts')

 load('generate_3d.mat','dU_ts_new','dE_ts_new','dN_ts_new');
   % Create selected data
   for c=1:length(index)
       dU_tes(:,c)=dU_ts_new(index(c,1),:);
       dE_tes(:,c)=dE_ts_new(index(c,1),:);
       dN_tes(:,c)=dN_ts_new(index(c,1),:);
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
   calc_std_3d=[mean_dU;std_dU;mean_dE;std_dE;mean_dN;std_dN];
   if exist('data.mat','file')
          save('calc_std.mat','calc_std_3d','-append');
   else
          save('calc_std.mat','calc_std_3d');
   end
   dlmwrite('3d_std.txt',calc_std_3d,'precision',8,'delimiter',' ');
   dlmwrite('index.txt',length(index));

   % Plot the mean scatter with error bar
   figure('rend','painters','pos',[10 10 800 600])
   
   % Plot dU
   subplot(3,1,1);
   plot(range,mean_dU,'k','LineWidth',1);
   hold on;
   errorbar(range,mean_dU,std_dU,'s','MarkerSize',1,'MarkerEdgeColor','red','MarkerFaceColor','red')
   datetick('x','mm-yy','keepticks')
   ylabel('dU disp (mm)')
   title('Error bar plot of 3D displacement result')
   
   % Plot dE
   subplot(3,1,2);
   plot(range,mean_dE,'k','LineWidth',1);
   hold on;
   errorbar(range,mean_dE,std_dE,'s','MarkerSize',1,'MarkerEdgeColor','red','MarkerFaceColor','red')
   datetick('x','mm-yy','keepticks')
   ylabel('dE disp (mm)')

   % Plot dN
   subplot(3,1,3);
   plot(range,mean_dN,'k','LineWidth',1);
   hold on;
   errorbar(range,mean_dN,std_dN,'s','MarkerSize',1,'MarkerEdgeColor','red','MarkerFaceColor','red')
   datetick('x','mm-yy','keepticks')
   ylabel('dN disp (mm)')
   xlabel('time')

   print(gcf,'3d_errorbar.png','-dpng','-r600');

end




