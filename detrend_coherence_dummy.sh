# 09.07.2018	NI	Calculate coherence on each master-slave pair
#			Calculate the mean coherence of SB network configuration
#			Take a threshold value of the mean coherence

# SET PARAMETER
list=$1
az_lks=1
dec_rng=1
wavelength=200
filter1=$GMTSAR/gmtsar/filters/gauss5x5
filter3=$GMTSAR/gmtsar/filters/fill.3x3
topo=$2 #/home/isya/APPS/ciloto/Sentinel1/batch_dsc/topo
###########################################
thresh=5.e-21
dec=1
scale=-JX6.5i

echo ""
echo "== Compute the Mean Coherence =="
shopt -s extglob
IFS=" "
num=1
while read master slave
do
  #cp ../raw_orig/$master.PRM .
  #cp ../raw_orig/$slave.PRM .

# make the custom filter2 and set the decimation
  make_gaussian_filter $master.PRM $dec_rng $az_lks $wavelength > ijdec
  filter2=gauss_$wavelength
  idec=$(cat ijdec | awk -v dc="$dec" '{ print dc*$1 }')
  jdec=$(cat ijdec | awk -v dc="$dec" '{ print dc*$2 }')
  echo "$filter2 $idec $jdec ($az_lks $dec_rng)"
  echo "idec = $idec , jdec = $jdec" 

  echo "making amplitudes..."
  conv $az_lks $dec_rng $filter1 $master.PRM amp1_tmp.grd=bf
  conv $idec $jdec $filter2 amp1_tmp.grd=bf amp1.grd
  rm amp1_tmp.grd
  conv $az_lks $dec_rng $filter1 $slave.PRM amp2_tmp.grd=bf
  conv $idec $jdec $filter2 amp2_tmp.grd=bf amp2.grd
  rm amp2_tmp.grd

  echo "making interferogram and filtering"
  intf.csh $master.PRM $slave.PRM -topo $topo/topo_ra.grd
  conv $az_lks $dec_rng $filter1 real.grd=bf real_tmp.grd=bf
  conv $idec $jdec $filter2 real_tmp.grd=bf realfilt.grd
  rm real_tmp.grd 
  rm real.grd
  conv $az_lks $dec_rng $filter1 imag.grd=bf imag_tmp.grd=bf
  conv $idec $jdec $filter2 imag_tmp.grd=bf imagfilt.grd
  rm imag_tmp.grd 
  rm imag.grd
  gmt grdmath realfilt.grd imagfilt.grd HYPOT  = amp.grd
  rm realfilt.grd imagfilt.grd

  echo "making correlation..."
  gmt grdmath amp1.grd amp2.grd MUL = tmp.grd
  gmt grdmath tmp.grd $thresh GE 0 NAN = mask.grd
  gmt grdmath amp.grd tmp.grd SQRT DIV mask.grd MUL FLIPUD = tmp2.grd=bf
  conv 1 1 $filter3 tmp2.grd=bf corr.grd

# calculate the mean coherence from each master-slave pair	
  if [ $num == 1 ]; then
    gmt grdmath corr.grd = sum.grd  
  else
    gmt grdmath sum.grd corr.grd ADD = sumtmp.grd 
    mv sumtmp.grd sum.grd 
  fi
  (( num++ ))

  rm amp1.grd amp2.grd amp.grd corr.grd
  #rm $master.PRM $slave.PRM $master.PRM0 $slave.PRM0
done < $list

  (( num-- )) 
  gmt grdmath sum.grd $num DIV = mean_corr.grd

# plot the average coherence
  gmt makecpt -T0./1/0.1 -Cgray -Z -N > corr.cpt
  echo "N  255   255   254" >> corr.cpt
  gmt grdimage mean_corr.grd $scale -Ccorr.cpt -Bf1000a5000WSen -X1.3 -Y3 -P -K > mean_corr.ps
  gmt psscale -Dx1.3/-1.5+w5/0.2+h+e -Ccorr.cpt -B0.2:correlation: -O >> mean_corr.ps
