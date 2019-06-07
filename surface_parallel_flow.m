function [ts]=surface_parallel_flow(value_type,varargin)

% 12.02.2019	NI	; Generate pseudo 3D vectors from LOS ascending, descending and
%                         surface-parallel-flow assumption
% Source: Joughin et al, 1998. Interferometric Estimation of Three-Dimensional Ice-Flow Using Ascending and Descending Passes, IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. 36.			  

% TYPE:
% surface_parallel_flow('ts') to calculate dU,dE,dN for time series
% surface_parallel_flow('mean') to calculate dU,dE,dN for the mean velocity

stdargin = nargin ;

if strcmp(value_type,'ts')

 X = sprintf('Calculate 3d vectors for time series using surface-parallel-flow + aspect');
   disp(X)

   %cd timeseries
   %load('data_match.mat','var_vector_asc','var_vector_dsc','var_lonlat')
   load('data_match.mat')
   if exist('uw_correct_asc')
      clear var_vector_asc var_vector_dsc;
      var_vector_asc=uw_correct_asc;
      var_vector_dsc=uw_correct_dsc;
      clear uw_correct_asc uw_correct_dsc;
   end
   pseudo=zeros(length(var_vector_asc),1);
   %load('../data_match.mat', 'var_angle', 'theta')

   load('interpolate.mat', 'range')
   dU_ts=zeros(length(var_angle),length(range));
   dE_ts=zeros(length(var_angle),length(range));
   dN_ts=zeros(length(var_angle),length(range));

   %% calculate dU,dE,dN for time series
	for n=1:length(range)
	    for c=1:length(var_vector_asc)
		% estimate dU,dE,dN with original least square
	    	A_mat=[var_vector_asc(c,n);var_vector_dsc(c,n);pseudo(c,1)];
	    	B1=cosd(var_angle(c,3));
		B2=-sind(var_angle(c,3)).*sind(var_angle(c,1)+90); 
		B3=-sind(var_angle(c,3)).*cosd(var_angle(c,1)+90);
		B4=cosd(var_angle(c,4));
		B5=-sind(var_angle(c,4)).*sind(var_angle(c,2)+90);
		B6=-sind(var_angle(c,4)).*cosd(var_angle(c,2)+90);

       		% based on surface flow direction
        	% define A
        	beta_angle(c,1)=-1*(var_angle(c,1));
        	alpha_angle(c,1)=180-((180+var_angle(c,2))+beta_angle(c,1));
        	A1=cosd(beta_angle(c,1));
        	A2=cosd(beta_angle(c,1)+alpha_angle(c,1));
        	A3=sind(beta_angle(c,1));
        	A4=sind(beta_angle(c,1)+alpha_angle(c,1));
        	A=[A1 A2;A3 A4];

        	% define B 
        	Be1=1/(sind(alpha_angle(c,1)).^2);
        	Be=Be1*[1 -cosd(alpha_angle(c,1));-cosd(alpha_angle(c,1)) 1];

        	% define C
        	C=[cotd(var_angle(c,3)) cotd(var_angle(c,3));cotd(var_angle(c,4)) cotd(var_angle(c,4))];
%        C=[-csc(var_angle(c,3)).^2 -csc(var_angle(c,3)).^2;-csc(var_angle(c,3)).^2 -csc(var_angle(c,3)).^2];

       		% define D
%               time_duration=range(n+1,1)-range(n,1); % in days
                time_duration=1;
        	% in radian
%        lambda=0.055465800000000; %Sentinel-1 wavelength in m
%        var_rad_1=var_vector(c,1)*-4*pi/lambda/1000;
%        var_rad_2=var_vector(c,2)*-4*pi/lambda/1000;
%        Da=var_rad_1/(time_asc*(4*pi/(lambda*1000))*(sind(var_angle(c,3))));
%        Dd=var_rad_2/(time_dsc*(4*pi/(lambda*1000))*(sind(var_angle(c,4))));
                % in dLOS
                Da=var_vector_asc(c,n)/(time_duration*(sind(var_angle(c,3))));
                Dd=var_vector_dsc(c,n)/(time_duration*(sind(var_angle(c,4))));
        	D=[Da;Dd];

        	% calculate vh in coordinates xy (H)
        	H=(1./(1-(A*Be*C)))*A*Be*D;

        	B7=-1;
        	B8=H(1,1);
        	B9=H(2,1);

		B=[B1 B2 B3;B4 B5 B6;B7 B8 B9];
	   	% calculate m --> [dU;dE;dN] vectors
	   	m(:,c)=lscov(B,A_mat);
	    end
	    dU_ts(:,n)=[m(1,:)'];
	    dE_ts(:,n)=[m(2,:)'];
	    dN_ts(:,n)=[m(3,:)'];
	    clear B1 B2 B3 B4 B5 B6 B7 B8 B9 c A B m;
	end
	clear n;

%        for n=1:length(range)
%            for c=1:length(var_angle)
%                if (var_angle(c,5) >= 0) && (var_angle(c,5) <= 180) && (dE_ts(c,n)<= 0)
%                   dN_tmp(c,n)=dN_ts(c,n)*-1;
%                   dN_ts(c,n)=dN_tmp(c,n);
%                elseif (var_angle(c,5) > 180) && (var_angle(c,5) <= 360) && (dE_ts(c,n)>= 0)
%                   dN_tmp(c,n)=dN_ts(c,n)*-1;
%                   dN_ts(c,n)=dN_tmp(c,n);
%                end
%            end
%        end

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

        % arrange dN based on aspect kuadrant to Y axis (from dN)
        for n=1:length(range)
            for c=1:length(var_angle)
                if (var_angle(c,5) >= 0) && (var_angle(c,5) <= 90) && (dN_ts_new(c,n)<= 0)
                   dN_tmp(c,n)=dN_ts_new(c,n)*-1;
                   dN_ts_new(c,n)=dN_tmp(c,n);
                elseif (var_angle(c,5) > 270) && (var_angle(c,5) <= 360) && (dN_ts_new(c,n)<= 0)
                   dN_tmp(c,n)=dN_ts_new(c,n)*-1;
                   dN_ts_new(c,n)=dN_tmp(c,n);
                elseif (var_angle(c,5) > 90) && (var_angle(c,5) <= 270) && (dN_ts_new(c,n)>= 0)
                   dN_tmp(c,n)=dN_ts_new(c,n)*-1;
                   dN_ts_new(c,n)=dN_tmp(c,n);
                end
            end
        end

        % see the vertical scale for plotting
        Umin=min(dU_ts_new);
        Umax=max(dU_ts_new);
        scale=[min(Umin) max(Umax)];

	date=datetime(range,'ConvertFrom','datenum');
	dates=datestr(date);
	dlmwrite('date.in',dates);
	clear date dates

	dlmwrite('dU_ts_3d_SF.txt',dU_ts_new,'precision',8,'delimiter',' ');
	dlmwrite('dE_ts_3d_SF.txt',dE_ts_new,'precision',8,'delimiter',' ');
	dlmwrite('dN_ts_3d_SF.txt',dN_ts_new,'precision',8,'delimiter',' ');
	dlmwrite('lonlat.txt',var_lonlat,'precision',8,'delimiter',' ');
        dlmwrite('ver_scale.txt',scale,'precision',8,'delimiter',' ');
	if exist('generate_3d_SF.mat','file')
	   save('generate_3d_SF.mat','dU_ts','-append');
	else
	   save('generate_3d_SF.mat','dU_ts');
	end
	save('generate_3d_SF.mat','dU_ts_new','-append');
	save('generate_3d_SF.mat','dE_ts','-append');
	save('generate_3d_SF.mat','dE_ts_new','-append');
	save('generate_3d_SF.mat','dN_ts','-append');
	save('generate_3d_SF.mat','dN_ts_new','-append');

   %cd ..

else

   X = sprintf('Calculate 3d vectors for the mean velocity (mm/year) with surface-parallel-flow');
   disp(X)

   load('data_match.mat')
   pseudo=zeros(length(var_vector),1);

%%%%%%%%% generate 3D %%%%%%%%%%%%%%%

   for c=1:length(var_vector)
    	% estimate dU,dE,dN with original least square (OLS)
    	A_mat=[var_vector(c,1);var_vector(c,2);pseudo(c,1)];
    	B1=cosd(var_angle(c,3));
	B2=-sind(var_angle(c,3)).*sind(var_angle(c,1)+90); 
	B3=-sind(var_angle(c,3)).*cosd(var_angle(c,1)+90);
	B4=cosd(var_angle(c,4));
	B5=-sind(var_angle(c,4)).*sind(var_angle(c,2)+90);
	B6=-sind(var_angle(c,4)).*cosd(var_angle(c,2)+90);

        % based on surface flow direction
        % define A
        beta_angle(c,1)=-1*(var_angle(c,1));
        alpha_angle(c,1)=180-((180+var_angle(c,2))+beta_angle(c,1));
        A1=cosd(beta_angle(c,1));
        A2=cosd(beta_angle(c,1)+alpha_angle(c,1));
        A3=sind(beta_angle(c,1));
        A4=sind(beta_angle(c,1)+alpha_angle(c,1));
        A=[A1 A2;A3 A4];

        % define B 
        Be1=1/(sind(alpha_angle(c,1)).^2);
        Be=Be1*[1 -cosd(alpha_angle(c,1));-cosd(alpha_angle(c,1)) 1];

        % define C
        C=[cotd(var_angle(c,3)) cotd(var_angle(c,3));cotd(var_angle(c,4)) cotd(var_angle(c,4))];
%        C=[-csc(var_angle(c,3)).^2 -csc(var_angle(c,3)).^2;-csc(var_angle(c,3)).^2 -csc(var_angle(c,3)).^2];

        % define D
        load('ps2_asc.mat','day')
        time_asc=(day(length(day),1)-day(1,1))/365.25; % in year
        clear day;
        load('ps2_dsc.mat','day')
        time_dsc=(day(length(day),1)-day(1,1))/365.25; % in year
        clear day;
        % in radian
%        lambda=0.055465800000000; %Sentinel-1 wavelength in m
%        var_rad_1=var_vector(c,1)*-4*pi/lambda/1000;
%        var_rad_2=var_vector(c,2)*-4*pi/lambda/1000;
%        Da=var_rad_1/(time_asc*(4*pi/(lambda*1000))*(sind(var_angle(c,3))));
%        Dd=var_rad_2/(time_dsc*(4*pi/(lambda*1000))*(sind(var_angle(c,4))));
        % in dLOS
%        Da=var_vector(c,1)/(time_asc*(sind(var_angle(c,3)))); --> if you want to calculate velocity !
%        Dd=var_vector(c,2)/(time_dsc*(sind(var_angle(c,4)))); --> the data is already in vel unit
        Da=var_vector(c,1)/(sind(var_angle(c,3)));
        Dd=var_vector(c,2)/(sind(var_angle(c,4)));
        D=[Da;Dd];

        % calculate vh in coordinates xy (H)
        H=(1./(1-(A*Be*C)))*A*Be*D;

        B7=-1;
        B8=H(1,1);
        B9=H(2,1);
        B=[B1 B2 B3;B4 B5 B6;B7 B8 B9];
	% calculate m --> [dU;dE;dN] vectors
	m(:,c)=lscov(B,A_mat);
   end

   dU=[var_lonlat(:,1) var_lonlat(:,2) m(1,:)'];
   dE=[var_lonlat(:,1) var_lonlat(:,2) m(2,:)'];
   dN=[var_lonlat(:,1) var_lonlat(:,2) m(3,:)'];

   clear B1 B2 B3 B4 B5 B6 B7 B8 B9 c A_mat B m;
   % see the vertical scale for plotting
   scale=[min(dU(:,3)) max(dU(:,3))];
   % save data
   dlmwrite('dE_3d_SF.txt',dE,'precision',8,'delimiter',' ');
   dlmwrite('dU_3d_SF.txt',dU,'precision',8,'delimiter',' ');
   dlmwrite('dN_3d_SF.txt',dN,'precision',8,'delimiter',' ');
   dlmwrite('ver_scale_SF.txt',scale,'precision',8,'delimiter',' ');
   if exist('generate_3d_SF.mat','file')
	   save('generate_3d_SF.mat','dU','-append');
   else
	   save('generate_3d_SF.mat','dU');
   end
   save('generate_3d_SF.mat','dE','-append');
   save('generate_3d_SF.mat','dN','-append');

end
