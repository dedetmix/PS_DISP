#!/bin/bash

### SET PARAMETER ###
region=106.9932/106.9945/-6.7075/-6.7055
path_ts=timeseries
#####################

mkdir -p $path_ts/velocity

# create CPT file for vertical displacement
min=$(awk '{print $1}' $path_ts/ver_scale.txt)
max=$(awk '{print $2}' $path_ts/ver_scale.txt)
gmt makecpt -T$min/$max/0.1 -N -I-Crainbow  > ver.cpt

for i in {1..88}
    do
	number=$(printf "%03d" "$i")
	rm -f dEN_velo_$number.txt dU_$number.txt dE_$number.txt dN_$number.txt

	## dU
	awk -v var="$i" '{ print $var }' $path_ts/dU_ts2.txt > tmp_dU.txt
	
	while read -r -u3 lon lat && read -r -u4 ver
	do
	  echo $lon$','$lat$','$ver >> dU_$number.txt
	done 3< $path_ts/lonlat2.txt 4<tmp_dU.txt

	gmt surface dU_$number.txt -R$region -I0.0000138889 -GdU_surface.grd -T0.25 -C0.1
	#gmt xyz2grd dU_$number.txt -Dm/m/m/1/0 -GdU.grd -I0.000046296 -R$region

	## dE and dN
	awk -v var="$i" '{ print $var }' $path_ts/dE_ts2.txt > tmp_dE.txt
	awk -v var="$i" '{ print $var }' $path_ts/dN_ts2.txt > tmp_dN.txt
	
	while read -r -u3 lon lat && read -r -u4 horE
	do
	  echo $lon$','$lat$','$horE >> dE_$number.txt
	done 3< $path_ts/lonlat2.txt 4<tmp_dE.txt

	while read -r -u3 lon lat && read -r -u4 horN
	do
	  echo $lon$','$lat$','$horN >> dN_$number.txt
	done 3< $path_ts/lonlat2.txt 4<tmp_dN.txt

	IFS=','
	while read -r -u3 e1 e2 e3 && read -r -u4 n1 n2 n3 
	do
	   echo $e1$'\t'$e2$'\t'$e3$'\t'$n3$'\t'0$'\t'0$'\t'0 >> dEN_velo_$number.txt 
	done 3<dE_$number.txt 4<dN_$number.txt

	rm -f tmp_dE.txt tmp_dN.txt tmp_dU.txt dU_$number.txt dE_$number.txt dN_$number.txt
	rm -f puncakpass_$number.ps

	gmt psxy Puncak_Pass_Landslide.xy -JM25 -B0.001 -P -Y3 -R$region -Wthinner,red -K > puncakpass_$number.ps
	gmt psclip Puncak_Pass_Landslide.xy -JM25 -B0.001 -P -R$region -O -K >> puncakpass_$number.ps
	gmt grdimage dU_surface.grd -J -Cver.cpt -B -P -R$region -O -K  >> puncakpass_$number.ps
	gmt psvelo dEN_velo_$number.txt -R -J -W0.5p,black -Se0.5/0/6 -O -K >> puncakpass_$number.ps
	gmt psclip -C -O -K >> puncakpass_$number.ps
	gmt psscale -DjCB+w4i/1c+o0/2c+h -Cver.cpt -R -J -Bx2+l"vertical disp (mm)" -O >> puncakpass_$number.ps
	gmt psconvert puncakpass_$number.ps -Tg -E526
	mv puncakpass_$number.png $path_ts/velocity/.

	rm -f dU_surface.grd dU_$number.grd dE_$number.grd dN_$number.grd dEN_velo_$number.txt puncakpass_$number.ps
	((i++)) 
done

convert -resize 50% -delay 20 -loop 0 $path_ts/velocity/puncakpass_{001..020}.png puncakpass_vel.gif

