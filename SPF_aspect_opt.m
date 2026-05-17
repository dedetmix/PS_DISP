function [ts]=SPF_aspect_opt(value_type,varargin)

% 12.02.2019	NI	; Generate pseudo 3D vectors from LOS ascending, descending and
%                         surface-parallel-flow assumption
% 05.07.2019	NI	; A control of SPF using slope aspect
% 15.07.2019	NI	; An optimal control of SPF using slope aspect and VCE algorithm
% 24.10.2019	NI	; Use with/without sensitivity
% Source: 
% Joughin et al, 1998. Interferometric Estimation of Three-Dimensional Ice-Flow Using Ascending and Descending Passes, IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. 36.
% Isya, N. H., Niemeier, W., and Gerke, M.: 3D ESTIMATION OF SLOW GROUND MOTION USING INSAR AND THE SLOPE ASPECT ASSUMPTION, A CASE STUDY: THE PUNCAK PASS LANDSLIDE, INDONESIA, ISPRS Ann. Photogramm. Remote Sens. Spatial Inf. Sci., IV-2/W5, 623-630, 2019.
%Jun, H.; Li, Z.-W.; Li, J.; Zhang, L.; Ding, X.; Jian-Jun, Z. & Sun, Q. 3-D movement mapping of the alpine glacier in Qinghai-Tibetan Plateau by integrating D-InSAR, MAI and Offset-Tracking: Case study of the Dongkemadi Glacier Global and Planetary Change, 2014, 118 			  

% TYPE:
% SPF_aspect_opt('ts') to calculate dU,dE,dN for time series
% SPF_aspect_opt('mean') to calculate dU,dE,dN for the mean velocity

stdargin = nargin ;

% define the time duration for SPF
fileID = fopen('spf_parm.txt');
spf_parm = textscan(fileID,'%s');
spf_parm = char(spf_parm{1});
fclose(fileID);

% define the sensitivity
fileID = fopen('vce_parm.txt');
vce_parm = textscan(fileID,'%s');
vce_parm = char(vce_parm{1});
fclose(fileID);

% define the convergent value
value_parm = dlmread('value_parm.txt');
value_sensitivity = value_parm(1);
value_convergent = value_parm(2);

if strcmp(value_type,'ts')
   
   %diary log_vce_ts.txt
   X = sprintf('Calculate 3d vectors for time series using surface-parallel-flow and aspect (optimal)');
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
   pseudo=zeros(length(theta),1);
   %load('../data_match.mat', 'var_angle', 'theta')

   load('interpolate.mat', 'range')
   dU_ts=zeros(length(var_angle),length(range));
   dE_ts=zeros(length(var_angle),length(range));
   dN_ts=zeros(length(var_angle),length(range));

   % adjust the first acq. time to be "0" value vor var_vector
   var_vector_asc_tmp=zeros(size(var_vector_asc));
   var_vector_dsc_tmp=zeros(size(var_vector_dsc));
   for n=1:length(range)
	    var_vector_asc_tmp(:,n)=var_vector_asc(:,n) - var_vector_asc(:,1);
	    var_vector_dsc_tmp(:,n)=var_vector_dsc(:,n) - var_vector_dsc(:,1);
   end
   clear n var_vector_asc var_vector_dsc;
   var_vector_asc=var_vector_asc_tmp;
   var_vector_dsc=var_vector_dsc_tmp;
   clear var_vector_asc_tmp var_vector_dsc_tmp;

   % The first inversion using a known weigthed matrix
   % Avoid P has a negative value for the first time series
   % Check R Index ascending
   for c=1:length(var_angle)
       if (rindex_asc(c) < 0)
          rindex_asc_tmp(c,1)=rindex_asc(c)*-1;
       else
          rindex_asc_tmp(c,1)=rindex_asc(c);
       end
   end
   rindex_asc=rindex_asc_tmp;
   clear rindex_asc_tmp 
   % Check R Index desscending
   for c=1:length(var_angle)
       if (rindex_dsc(c) < 0)
          rindex_dsc_tmp(c,1)=rindex_dsc(c)*-1;
       else
          rindex_dsc_tmp(c,1)=rindex_dsc(c);
       end
   end
   rindex_dsc=rindex_dsc_tmp;
   clear rindex_dsc_tmp

   diff_r=rindex_asc-rindex_dsc;

   %% calculate dU,dE,dN for time series
	for n=1:length(range)
	    for c=1:length(theta)
                    if (rem(c,500) == 0)
                       A = sprintf('Date: %i for Pixel : %i',n,c);
                       disp(A)
                    end
		% estimate dU,dE,dN with original least square
	    	A_mat=[var_vector_asc(c,n);var_vector_dsc(c,n);pseudo(c,1);pseudo(c,1)];
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

                % define duration based on linear/non-linear assumption
   		if strcmp(spf_parm,'non-linear')
         	   time_duration=1;
   		else
      		   if (n == 1)
         	      time_duration=1;
      		   else
         	      time_duration=range(n)-range(1); % in days
      		   end
   		end

   		% define with/without sensitivity for D and VCE
   		if strcmp(vce_parm,'sensitivity')
      		   P=[rindex_asc(c) rindex_dsc(c) 0.41 0.155];
   		else
      		   P=[1 1 1 1];
   		end

      		m3=0.41*length(range);
      		m4=0.115*length(range);

                % define D   
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

                % based on slope aspect
		B10=0;
		B11=cosd(theta(c,1))*cosd(90-theta(c,1)); %scalar projection from dE
		B12=-1;

		B=[B1 B2 B3;B4 B5 B6;B7 B8 B9;B10 B11 B12];

	   	% calculate m --> [dU;dE;dN] vectors with VCE and iteration
                GLS_error=1;
                iter=1;
                while GLS_error >= value_convergent
                  vce;
                  % save statistic result
                  vce_sta{c,n}=computation;
                  if (rem(c,500) == 0)
                     X = sprintf('%i iteration to estimate the weight matrices',iter);
                     disp(X)
                     Y = sprintf('GLS error = %f',GLS_error);
                     disp(Y)
                  end
                  iter=iter+1;
                  m(:,c)=m_tmp;
                end
                clear computation;

	    end
	    dU_ts(:,n)=[m(1,:)'];
	    dE_ts(:,n)=[m(2,:)'];
	    dN_ts(:,n)=[m(3,:)'];
	    clear B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 c A B m;
	end
	clear n;

       % arrange dN based on aspect kuadrant to Y axis (from dE)
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
%	dU_ts_new=zeros(size(dU_ts));
%	dE_ts_new=zeros(size(dE_ts));
%	dN_ts_new=zeros(size(dN_ts));
%	for n=1:length(range)
%	    dU_ts_new(:,n)=dU_ts(:,n) - dU_ts(:,1);
%	    dE_ts_new(:,n)=dE_ts(:,n) - dE_ts(:,1);
%	    dN_ts_new(:,n)=dN_ts(:,n) - dN_ts(:,1);
%	end
%	clear n;

        % save the original scalar projection of dN
        dU_ts_new=dU_ts;
        dE_ts_new=dE_ts;
        dN_ts_new=dN_ts;

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

	dlmwrite('dU_ts_3d_SF_opt.txt',dU_ts_new,'precision',8,'delimiter',' ');
	dlmwrite('dE_ts_3d_SF_opt.txt',dE_ts_new,'precision',8,'delimiter',' ');
	dlmwrite('dN_ts_3d_SF_opt.txt',dN_ts_new,'precision',8,'delimiter',' ');
	dlmwrite('lonlat.txt',var_lonlat,'precision',8,'delimiter',' ');
        dlmwrite('ver_scale_3d_SF_opt.txt',scale,'precision',8,'delimiter',' ');
	if exist('generate_3d_SF_opt.mat','file')
	   save('generate_3d_SF_opt.mat','dU_ts','-append');
	else
	   save('generate_3d_SF_opt.mat','dU_ts');
	end
	save('generate_3d_SF_opt.mat','dU_ts_new','-append');
	save('generate_3d_SF_opt.mat','dE_ts','-append');
	save('generate_3d_SF_opt.mat','dE_ts_new','-append');
	save('generate_3d_SF_opt.mat','dN_ts','-append');
	save('generate_3d_SF_opt.mat','dN_ts_new','-append');
        % save statistical computation
        colNames = {'GLS_error','P_grup_1','P_grup_2','P_grup_3','P_grup_4','rmse','std_error_dU','std_error_dE','std_error_dN'};
        rowNames = {'iteration'};
	%if exist('vce_SF_opt.mat','file')
	%   save('vce_SF_opt.mat','vce_sta','-append');
	%else
	   save('vce_SF_opt.mat','vce_sta');
	%end
	save('vce_SF_opt.mat','colNames','-append');
	save('vce_SF_opt.mat','rowNames','-append');
   
   %diary off
   %cd ..

else

   X = sprintf('Step is unknown');
   disp(X)

end
