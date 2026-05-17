#!/bin/bash

mode=$1

dir=$(pwd)
path=$dir #asc_surface_V-DAS_v2 #dir
psfile=dEN
psfile_gps=gps_hor.ps
path_ts=timeseries
region=106.99/107.02/-6.73/-6.7

mkdir -p $path_ts/horizontal
mkdir -p $path_ts/horizontal/dEN_velo

if [ "$mode" == "ts" ]; then #use "$variable" if defining string format

   echo "== display time series (horizontal) ==="
   echo " "

   for i in {1..88}
   do
	number=$(printf "%03d" "$i")
	awk -v var="$i" '{ print $var }' $path_ts/dE_ts.txt > tmp_dE.txt
	awk -v var="$i" '{ print $var }' $path_ts/dN_ts.txt > tmp_dN.txt
	
	while read -r -u3 lon lat && read -r -u4 horE
	do
	  echo $lon$','$lat$','$horE >> dE_$number.txt
	done 3< $path_ts/lonlat.txt 4<tmp_dE.txt

	while read -r -u3 lon lat && read -r -u4 horN
	do
	  echo $lon$','$lat$','$horN >> dN_$number.txt
	done 3< $path_ts/lonlat2.txt 4<tmp_dN.txt

	rm -f dEN_velo_$number.txt
	gmt xyz2grd dE_$number.txt -Dm/m/m/1/0 -GdE.grd -I0.000046296 -R$region
	gmt xyz2grd dN_$number.txt -Dm/m/m/1/0 -GdN.grd -I0.000046296 -R$region
	gmt grdsample dE.grd -I0.000138889 -GdE_re.grd
	gmt grdsample dN.grd -I0.000138889 -Gtmp.grd
	gmt grdclip tmp.grd -GdN_re.grd -Sa10/NaN -Sb-10/NaN -V
	rm tmp.grd
	gmt grd2xyz -R$region dE_re.grd > dE_re.txt
	gmt grd2xyz -R$region dN_re.grd > dN_re.txt
	while read -r -u3 e1 e2 e3 && read -r -u4 n1 n2 n3 
	do
	  echo $e1$'\t'$e2$'\t'$e3$'\t'$n3$'\t'0$'\t'0$'\t'0 >> dEN_velo_$number.txt 
	done 3<dE_re.txt 4<dN_re.txt
 	rm -f dE_re.txt dE.grd dE_re.grd dE_$number.txt tmp_dE.txt
	rm -f dN_re.txt dN.grd dN_re.grd dN_$number.txt tmp_dN.txt

	date=$(sed "${i}q;d" $path_ts/date.in)
	gmt psxy ciloto_prone_landslide.xy -JX9i -R$region -K -Wthinner,red -B0.05 -Bwesn+t$date > $path/$psfile"_"$number.ps #--MAP_FRAME_TYPE=outside
	gmt psclip ciloto_prone_landslide.xy -J -R -O -K -V >> $path/$psfile"_"$number.ps
	gmt psvelo $path/dEN_velo_$number.txt -R -J -W0.5p,green -Se0.05/0/6 -V -O -K >> $path/$psfile"_"$number.ps
	gmt psclip -J -R -C -O -V >> $path/$psfile"_"$number.ps
	gmt psconvert $psfile"_"$number.ps -Tg -E526 -W+k+t"horizontal disp"+l256/-1 -V
	rm $path/$psfile"_"$number.ps
	mv dEN_velo_$number.txt $path_ts/horizontal/dEN_velo/.
	mv $path/$psfile"_"$number.png $path_ts/horizontal/.
	mv $path/$psfile"_"$number.kml $path_ts/horizontal/.

	((i++)) 
   done

   convert -resize 50% -delay 20 -loop 0 $path_ts/horizontal/dEN_{001..020}.png dEN.gif

else

	# create dEN_velo.txt
	rm -f dEN_velo.txt
	gmt xyz2grd dN.txt -Dm/m/m/1/0 -GdN.grd -I0.000046296 -R$region
	gmt grdsample dN.grd -I0.000138889 -GdN_re.grd
	gmt xyz2grd dE.txt -Dm/m/m/1/0 -GdE.grd -I0.000046296 -R$region
	gmt grdsample dE.grd -I0.000138889 -GdE_re.grd
	gmt grd2xyz -R$region dE_re.grd > dE_re.txt
	gmt grd2xyz -R$region dN_re.grd > dN_re.txt
	while read -r -u3 e1 e2 e3 && read -r -u4 n1 n2 n3 
	do
	echo $e1$'\t'$e2$'\t'$e3$'\t'$n3$'\t'0$'\t'0$'\t'0 >> dEN_velo.txt 
	done 3<dE_re.txt 4<dN_re.txt
	rm dE.grd dN.grd dN_re.grd dE_re.grd dE_re.txt dN_re.txt	

	# create map velocity for the horizontal vector
	gmt psxy ciloto_prone_landslide.xy -JX10i -R$region -K -Wthinner,red -B0.05wesn --MAP_FRAME_TYPE=inside > $path/$psfile.ps
	gmt psclip ciloto_prone_landslide.xy -J -R -O -K -V >> $path/$psfile.ps
	gmt psxy -J -R -B -W0.4p,purple -O -K profil_puncakpass.txt >> $path/$psfile.ps
	gmt psxy -J -R -B -Sc0.01i -Gpurple -O -K profil_puncakpass.txt >> $path/$psfile.ps
	gmt psvelo $path/dEN_velo.txt -R -J -W0.5p,orange -Se0.05/0/6 -V -O -K >> $path/$psfile.ps
	gmt psclip -J -R -C -O -V >> $path/$psfile.ps

	#convert ps to kml file
	cd $path
	gmt psconvert $psfile.ps -TG -W+k+t"horizontal disp"+l256/-1 -V -E526 -V #-Tg for white PNG, -TG for tranparent
	cd $dir

	##gmt psxy ciloto_prone_landslide.xy -JX10d -R -K -Wthinner,red -Bwesn > $psfile_gps
	##gmt psvelo gps_all_hor_geo.txt -R -J -W0.5p,black -Se0.25/0/6 -V -O >> $psfile_gps
	##gmt psconvert $psfile_gps -Tg -W+k+t"gps horizontal disp"+l256/-1 -V

fi
