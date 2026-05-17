% simulation SMALL area

% prepare R-Index from vslope calculation
load('../data.mat', 'los_asc')
load('../data.mat', 'los_dsc')  
load('../V_slope.mat', 'r_index_asc')
load('../V_slope.mat', 'r_index_dsc') 
       i=1;
       for c=1:length(los_asc)
           if (~isnan(los_asc(c,3))) && (~isnan(los_dsc(c,3)))
               rindex_asc(i,:)=[r_index_asc(c,1)];
               rindex_dsc(i,:)=[r_index_dsc(c,1)];
               i=i+1;
           end
       end
       clear c i;

% prepare ROI area
load('data_match.mat')
load('../data_match.mat','var_slope')

Blat=find(var_lonlat(:,2)>-6.714166 & var_lonlat(:,2)<-6.713736);
Blon=find(var_lonlat(:,1)>107.005786 & var_lonlat(:,1)<107.006292);
Blast=ismember(Blat,Blon);
indexes=find(Blast);
index=Blat(indexes);
for c=1:length(index)
    var_vector_asc_tmp(:,c)=var_vector_asc(index(c,1),:);
    var_vector_dsc_tmp(:,c)=var_vector_dsc(index(c,1),:);
    var_lonlat_tmp(:,c)=var_lonlat(index(c,1),:);
    var_angle_tmp(:,c)=var_angle(index(c,1),:);
    var_slope_tmp(:,c)=var_slope(index(c,1),:);
    rindex_asc_tmp(:,c)=rindex_asc(index(c,1),:);
    rindex_dsc_tmp(:,c)=rindex_dsc(index(c,1),:);
    theta_tmp(:,c)=theta(index(c,1),:);
end
var_vector_asc=var_vector_asc_tmp';
var_vector_dsc=var_vector_dsc_tmp';
var_lonlat=var_lonlat_tmp';
var_angle=var_angle_tmp';
var_slope=var_slope_tmp';
rindex_asc=rindex_asc_tmp';
rindex_dsc=rindex_dsc_tmp';
theta=theta_tmp';

% save data_match_tmp.mat
       if exist('data_match_tmp.mat','file')
          save('data_match_tmp.mat','var_angle','-append');
       else
          save('data_match_tmp.mat','var_angle');
       end
       save('data_match_tmp.mat','var_vector_asc','-append');
       save('data_match_tmp.mat','var_vector_dsc','-append');
       save('data_match_tmp.mat','var_lonlat','-append');
       save('data_match_tmp.mat','var_slope','-append');
       save('data_match_tmp.mat','rindex_asc','-append');
       save('data_match_tmp.mat','rindex_dsc','-append');
       save('data_match_tmp.mat','theta','-append');
