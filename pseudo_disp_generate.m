function [ts]=pseudo_disp_generate(value_type,varargin)

% 23.05.2018	NI	; Generate pseudo 3D vectors from LOS ascending, descending and aspect
%			  Previous parameters were generated from PS_DISP_matlab.m or disp_generate.m 

% TYPE:
% pseudo_disp_generate('ts') to calculate dU,dE,dN for time series
% pseudo_disp_generate('mean') to calculate dU,dE,dN for the mean velocity

stdargin = nargin ;

if strcmp(value_type,'ts')
   X = sprintf('Calculate 3d vectors for time series');
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
	    	A=[var_vector_asc(c,n);var_vector_dsc(c,n);pseudo(c,1)];
	    	B1=cosd(var_angle(c,3));
		B2=-sind(var_angle(c,3)).*sind(var_angle(c,1)+90); 
		B3=-sind(var_angle(c,3)).*cosd(var_angle(c,1)+90);
		B4=cosd(var_angle(c,4));
		B5=-sind(var_angle(c,4)).*sind(var_angle(c,2)+90);
		B6=-sind(var_angle(c,4)).*cosd(var_angle(c,2)+90);
		B7=0;
		  if (var_angle(c,5) >= 0) && (var_angle(c,5) <= 90)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
		  elseif (var_angle(c,5) > 90) && (var_angle(c,5) <= 180)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
		  elseif (var_angle(c,5) > 180) && (var_angle(c,5) <= 270)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
		  elseif (var_angle(c,5) > 270) && (var_angle(c,5) <= 360)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
		  elseif (var_angle(c,5) > 360)
			X = sprintf('%i theta set to quadrant I',c);
			disp(X)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
		  else
			X = sprintf('%i theta set to quadrant IV',c);
			disp(X)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
		  end
		%B8=cosd(theta(c,1))*cosd(90-theta(c,1)); %see disp_generate.m to preview "theta" 
		B9=-1;
		B=[B1 B2 B3;B4 B5 B6;B7 B8 B9];
	   	% calculate m --> [dU;dE;dN] vectors
	   	m(:,c)=lscov(B,A);
	    end
	    dU_ts(:,n)=[m(1,:)'];
	    dE_ts(:,n)=[m(2,:)'];
	    dN_ts(:,n)=[m(3,:)'];
	    clear B1 B2 B3 B4 B5 B6 B7 B8 B9 c A B m;
	end
	clear n;

        % arrange dN based on aspect kuadrant
        for n=1:length(range)
            for c=1:length(var_angle)
                if (var_angle(c,5) >= 0) && (var_angle(c,5) <= 180) && (dE_ts(c,n)<= 0)
                   dN_tmp(c,n)=dN_ts(c,n)*-1;
                   dN_ts(c,n)=dN_tmp(c,n);
                elseif (var_angle(c,5) >= 180) && (var_angle(c,5) <= 360) && (dE_ts(c,n)>= 0)
                   dN_tmp(c,n)=dN_ts(c,n)*-1;
                   dN_ts(c,n)=dN_tmp(c,n);
                end
            end
         end

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

        % see the vertical scale for plotting
        Umin=min(dU_ts_new);
        Umax=max(dU_ts_new);
        scale=[min(Umin) max(Umax)];

	date=datetime(range,'ConvertFrom','datenum');
	dates=datestr(date);
	dlmwrite('date.in',dates);
	clear date dates

	dlmwrite('dU_ts.txt',dU_ts_new,'precision',8,'delimiter',' ');
	dlmwrite('dE_ts.txt',dE_ts_new,'precision',8,'delimiter',' ');
	dlmwrite('dN_ts.txt',dN_ts_new,'precision',8,'delimiter',' ');
	dlmwrite('lonlat.txt',var_lonlat,'precision',8,'delimiter',' ');
        dlmwrite('ver_scale.txt',scale,'precision',8,'delimiter',' ');
	if exist('generate_3d_pseudo.mat','file')
	   save('generate_3d_pseudo.mat','dU_ts','-append');
	else
	   save('generate_3d_pseudo.mat','dU_ts');
	end
	save('generate_3d_pseudo.mat','dU_ts_new','-append');
	save('generate_3d_pseudo.mat','dE_ts','-append');
	save('generate_3d_pseudo.mat','dE_ts_new','-append');
	save('generate_3d_pseudo.mat','dN_ts','-append');
	save('generate_3d_pseudo.mat','dN_ts_new','-append');

   cd ..
else
   
   X = sprintf('Calculate 3d vectors for the mean velocity (mm/year)');
   disp(X)

   load('data_match.mat')
   pseudo=zeros(length(var_vector),1);

   for c=1:length(var_vector)
    	% estimate dU,dE,dN with original least square (OLS)
    	A=[var_vector(c,1);var_vector(c,2);pseudo(c,1)];
    	B1=cosd(var_angle(c,3));
	B2=-sind(var_angle(c,3)).*sind(var_angle(c,1)+90); 
	B3=-sind(var_angle(c,3)).*cosd(var_angle(c,1)+90);
	B4=cosd(var_angle(c,4));
	B5=-sind(var_angle(c,4)).*sind(var_angle(c,2)+90);
	B6=-sind(var_angle(c,4)).*cosd(var_angle(c,2)+90);
	B7=0;
	if (var_angle(c,5) >= 0) && (var_angle(c,5) <= 90)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
	elseif (var_angle(c,5) > 90) && (var_angle(c,5) <= 180)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
	elseif (var_angle(c,5) > 180) && (var_angle(c,5) <= 270)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
	elseif (var_angle(c,5) > 270) && (var_angle(c,5) <= 360)
			B8=cosd(theta(c,1))*cosd(90-theta(c,1));
        elseif (var_angle(c,5) > 360)
	       X = sprintf('%i theta set to quadrant I',c);
	       disp(X)
	       B8=cosd(theta(c,1))*cosd(90-theta(c,1));
	else
	       X = sprintf('%i theta set to quadrant IV',c);
	       disp(X)
	       B8=cosd(theta(c,1))*cosd(90-theta(c,1));
	end
	%B8=cosd(theta(c,1))*cosd(90-theta(c,1)); %see disp_generate.m to preview "theta" 
	B9=-1;
	B=[B1 B2 B3;B4 B5 B6;B7 B8 B9];
	% calculate m --> [dU;dE;dN] vectors
	m(:,c)=lscov(B,A);
   end
       
   dU=[var_lonlat(:,1) var_lonlat(:,2) m(1,:)'];
   dE=[var_lonlat(:,1) var_lonlat(:,2) m(2,:)'];
   dN=[var_lonlat(:,1) var_lonlat(:,2) m(3,:)'];
   % arrange dN based on aspect kuadran
   for c=1:length(var_angle)
       if (var_angle(c,5) >= 0) && (var_angle(c,5) <= 180) && (dE(c,3)<= 0)
          dN_tmp(c,1)=dN(c,3)*-1;
          dN(c,3)=dN_tmp(c,1);
	  X = sprintf('%i fix the aspect direction for dN',c);
	  disp(X)
       elseif (var_angle(c,5) >= 180) && (var_angle(c,5) <= 360) && (dE(c,3)>= 0)
          dN_tmp(c,1)=dN(c,3)*-1;
          dN(c,3)=dN_tmp(c,1);
	  X = sprintf('%i fix the aspect direction for dN',c);
	  disp(X)
       end
   end
   clear B1 B2 B3 B4 B5 B6 B7 B8 B9 c A B m dN_tmp;

   % see the vertical scale for plotting
   scale=[min(dU(:,3)) max(dU(:,3))];
   % save data
   dlmwrite('dE.txt',dE,'precision',8,'delimiter',' ');
   dlmwrite('dU.txt',dU,'precision',8,'delimiter',' ');
   dlmwrite('dN.txt',dN,'precision',8,'delimiter',' ');
   dlmwrite('ver_scale.txt',scale,'precision',8,'delimiter',' ');
   if exist('generate_3d_pseudo.mat','file')
	   save('generate_3d_pseudo.mat','dU','-append');
   else
	   save('generate_3d_pseudo.mat','dU');
   end
   save('generate_3d_pseudo.mat','dE','-append');
   save('generate_3d_pseudo.mat','dN','-append');
end  
