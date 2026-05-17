% add aditional data for vce computation : R Index and Slope Angle

load('../V_slope.mat', 'r_index_asc')
load('../V_slope.mat', 'r_index_dsc')
load('../data_match.mat','var_slope')
load('../data.mat', 'az_angle_asc', 'az_angle_dsc', 'inc_angle_asc', 'inc_angle_dsc','aspect')
load('interpolate.mat', 'asc_nn')
load('interpolate.mat', 'dsc_nn')

i=1;
       for c=1:length(asc_nn)
           if (~isnan(asc_nn(c,1))) && (~isnan(dsc_nn(c,1))) && (~isnan(az_angle_asc(c,3))) && (~isnan(az_angle_dsc(c,3))) && (~isnan(inc_angle_asc(c,3))) && (~isnan(inc_angle_dsc(c,3))) && (~isnan(aspect(c,3)))
 rindex_asc(i,:)=[r_index_asc(c,1)];
               rindex_dsc(i,:)=[r_index_dsc(c,1)];
               i=i+1;
           end
       end
       clear c i;


save('data_match.mat','var_slope','-append');
save('data_match.mat','rindex_asc','-append');
save('data_match.mat','rindex_dsc','-append');
