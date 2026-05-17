% VCE algorithm

% 1) classify into 4 groups and open range for K option
     %P=[1 1 1 1];
     %load('interpolate.mat', 'range')
     %range=1;

% 2) do the first WLS solution
  [m_tmp,std_err,mse]=lscov(B,A_mat,P);

  %% calculate Generalized Survey Adjustment
  % define V and N
  V1=[B1 B2 B3]*m_tmp-A_mat;
  V2=[B4 B5 B6]*m_tmp-A_mat;
  V3=[B7 B8 B9]*m_tmp-A_mat;
  V4=[B10 B11 B12]*m_tmp-A_mat;
  V=V1+V2+V3+V4;
  N1=[B1 B2 B3]'*P(1)*[B1 B2 B3];
  N2=[B4 B5 B6]'*P(2)*[B4 B5 B6];
  N3=[B7 B8 B9]'*P(3)*[B7 B8 B9];
  N4=[B10 B11 B12]'*P(4)*[B10 B11 B12];
  N=N1+N2+N3+N4;
  % define W_the and S
  W_thet=[V1'*P(1)*V1; V2'*P(2)*V2; V3'*P(3)*V3; V4'*P(4)*V4];
  tmp_03=round(P(3),2);
  tmp_04=round(P(4),2);
  S_tmp1=length(range)-(2*trace(inv(N)*N1))+(trace(inv(N)*N1)); % m = the number of observations in the nth group (based on acq date)
  S_tmp2=length(range)-(2*trace(inv(N)*N2))+(trace(inv(N)*N2));
  S_tmp3=m3-(2*trace(inv(N)*N3))+(trace(inv(N)*N3)); % m = due to an assumption, taking from the average asc and dsc track
  S_tmp4=m4-(2*trace(inv(N)*N4))+(trace(inv(N)*N4));
  S=[S_tmp1^2 trace(inv(N)*N1*inv(N)*N2) trace(inv(N)*N1*inv(N)*N3) trace(inv(N)*N1*inv(N)*N4);%
     trace(inv(N)*N1*inv(N)*N2) S_tmp2^2 trace(inv(N)*N2*inv(N)*N3) trace(inv(N)*N2*inv(N)*N4);%
     trace(inv(N)*N1*inv(N)*N3) trace(inv(N)*N2*inv(N)*N3) S_tmp3^2 trace(inv(N)*N3*inv(N)*N4);%
     trace(inv(N)*N1*inv(N)*N4) trace(inv(N)*N2*inv(N)*N4) trace(inv(N)*N3*inv(N)*N4) S_tmp4^2];
  % calculate GLS_thet
  %GLS_thet=inv(S)*W_thet;
  GLS_thet=lscov(S,W_thet);
  % calculate the new P
  % define P based on R Index Asc or Dsc, the first P is based on low R value
  if (diff_r(c) <= 0.2)
     P1_tmp=P(1);
     P2_tmp=GLS_thet(1)/(GLS_thet(2)*inv(P(2)));
     P3_tmp=GLS_thet(1)/(GLS_thet(3)*inv(P(3)));
     P4_tmp=GLS_thet(1)/(GLS_thet(4)*inv(P(4)));
  else
     P1_tmp=P(2);
     P2_tmp=GLS_thet(2)/(GLS_thet(1)*inv(P(1)));
     P3_tmp=GLS_thet(2)/(GLS_thet(3)*inv(P(3)));
     P4_tmp=GLS_thet(2)/(GLS_thet(4)*inv(P(4)));
  end

% 3) replace the new P to be the first P and do iteration until the condition is fulfilled.
  P=[P1_tmp P2_tmp P3_tmp P4_tmp];
  GLS_error=max(GLS_thet)-min(GLS_thet);

% save statistic result
  computation(iter,:)=[GLS_error P(1) P(2) P(3) P(4) sqrt(mse) std_err(1) std_err(2) std_err(3)];
