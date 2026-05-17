% Optimal surface parallel flow

% Prepare Data
Blat=find(var_lonlat(:,2)>-6.714166 & var_lonlat(:,2)<-6.713736);
Blon=find(var_lonlat(:,1)>107.005786 & var_lonlat(:,1)<107.006292);
Blast=ismember(Blat,Blon);
indexes=find(Blast);
index=Blat(indexes);
for c=1:length(index)
    los_asc_tes(:,c)=var_vector_asc(index(c,1),:);
    los_dsc_tes(:,c)=var_vector_dsc(index(c,1),:);
    var_lonlat_tes(:,c)=var_lonlat(index(c,1),:);
    var_angle_tes(:,c)=var_angle(index(c,1),:);
end
los_asc_tes=los_asc_tes';
los_dsc_tes=los_dsc_tes';
var_lonlat_tes=var_lonlat_tes';
var_angle_tes=var_angle_tes';

%clear los_asc_tes los_dsc_tes var_lonlat_tes var_angle_tes Blat Blon Blast indexes index c 

% create convolution design matrix
l=(rows-3)+1;
c=(columns-3)+1;

for t_l=0:(l-1)
    n_l=t_l+1;
    for t_c=0:(c-1)
        n_c=t_c+1;
        tmp_matrix{n_l,n_c}=[asc(columns*(t_l)+(t_c+1)) asc(columns*(t_l)+(t_c+2)) asc(columns*(t_l)+(t_c+3)); asc(columns*(t_l+1)+(t_c+1)) asc(columns*(t_l+1)+(t_c+2)) asc(columns*(t_l+1)+(t_c+3)); asc(columns*(t_l+2)+(t_c+1)) asc(columns*(t_l+2)+(t_c+2)) asc(columns*(t_l+2)+(t_c+3))];
    end 
end

% arrange to be one vector column
tmp_Bcv=tmp_matrix{1,1}';
tmp_Bcv=tmp_Bcv(:);


