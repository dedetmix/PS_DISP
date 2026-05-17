#calculate aspect
gmt grdgradient dem.grd -Da -Gdem_aspect.grd
gmt grd2xyz dem_aspect.grd -R90.54/90.56/30.43/30.46 > dem_aspect.txt

#calculate slope
gmt grdgradient dem.grd -D -Gdem_slope.grd -Sslope.grd
grdmath slope.grd ATAN PI DIV 180 MUL = slope2.grd
mv slope2.grd dem_slope.grd
rm -f slope.grd

#export slope to matlab
# gmt xyz2grd dem_slope.txt -Ddegree/degree/degree/1/0 -Gslope.grd -I$resolution -R$region -V #optional if you have txt file of slope
gmt grdsample dem_slope.grd -I$resolution -Gslope_re.grd
gmt grd2xyz -R$region slope_re.grd > slope.xyz
# in Matlab
>> delimiterIn='\t';
>> slope=importdata('slope.xyz',delimiterIn);
>> save('data.mat','slope','-append');

# now work in matlab to use vslope.m to difine R index and C coeff
>> vslope('asc')
>> vslope('dsc')
>> load ('data.mat')
>> load('V_slope.mat', 'r_index_asc')
>> load('V_slope.mat', 'r_index_dsc')

       i=1;
       for c=1:length(los_asc)
           if (~isnan(los_asc(c,3))) && (~isnan(los_dsc(c,3))) && (~isnan(az_angle_asc(c,3))) && (~isnan(az_angle_dsc(c,3))) && (~isnan(inc_angle_asc(c,3))) && (~isnan(inc_angle_dsc(c,3))) && (~isnan(aspect(c,3))) 
              var_slope(i,:)=[slope(c,3)];
              rindex_asc(i,:)=[r_index_asc(c,1)];
              rindex_dsc(i,:)=[r_index_dsc(c,1)];
               i=i+1;
           end
       end
       clear c i;

save('data_match.mat','var_slope','-append');
save('data_match.mat','rindex_asc','-append');
save('data_match.mat','rindex_dsc','-append');
% on timeseries dir too
save('timeseries/data_match.mat','var_slope','-append');
save('timeseries/data_match.mat','rindex_asc','-append');
save('timeseries/data_match.mat','rindex_dsc','-append');
