#!/bin/bash

t=95
scale_length=80
method=SF_opt

gmt grdtrack -Gdem_re.grd lonlat.txt > lonlatelev.txt

i=1 #i=1 for time series

while [ $i -le $t ] 
do
   # extract dE
   awk -v var="$i" '{ print $var }' dE_ts_3d_"$method".txt > tmp_dE.txt
   # extract dN
   awk -v var="$i" '{ print $var }' dN_ts_3d_"$method".txt > tmp_dN.txt

   # calculate polar coordinates
   gmt math tmp_dE.txt tmp_dN.txt R2 = tmp1.txt
   gmt math tmp1.txt SQRT = tmp_dR_$i.txt
   gmt math tmp_dR_$i.txt $scale_length DIV = dR_$i.txt
   rm -f tmp1.txt tmp_dR_$i.txt
   gmt math tmp_dN.txt tmp_dE.txt ATAN2 = tmp_theta_$i.txt
   gmt math tmp_theta_$i.txt R2D = theta_$i.txt
   rm -f tmp_dE.txt tmp_dN.txt tmp_theta_$i.txt

   # arrange new file with 5 parameters
   # x	y	z	angle	length
     rm -f velo_hor_$i.txt
     IFS="\t"
     while read -r -u3 lon lat elev && read -r -u4 angle && read -r -u5 length
     do
       echo -e $lon"\t"$lat"\t"$elev"\t"$angle"\t"$length >> velo_hor_$i.txt
     done 3< lonlatelev.txt 4<theta_$i.txt 5<dR_$i.txt
   ((i++)) 
done
