#!/bin/bash

mode=$1

dir=$(pwd)
path=$dir #asc_surface_V-DAS_v2 or $dir
psfile=dU
psfile_gps=gps_ver.ps
path_ts=timeseries
region=106.99/107.02/-6.73/-6.7

mkdir -p $path_ts/vertical
mkdir -p $path_ts/vertical/dU_velo

if [ "$mode" == "ts" ]; then #use "$variable" if defining string format

   echo "== plot time series (vertical) ==="
   echo " "

   for i in {1..88}
   do
	number=$(printf "%03d" "$i")
	awk -v var="$i" '{ print $var }' $path_ts/dU_ts.txt > tmp_dU.txt
	
	while read -r -u3 lon lat && read -r -u4 ver
	do
	  echo $lon$','$lat$','$ver >> dU_$number.txt
	done 3< $path_ts/lonlat.txt 4<tmp_dU.txt

	rm -f dU_velo_$number.txt
	gmt xyz2grd dU_$number.txt -Dm/m/m/1/0 -GdU.grd -I0.000046296 -R$region
	gmt grdsample dU.grd -I0.000138889 -GdU_re.grd
	gmt grd2xyz -R$region dU_re.grd > dU_re.txt
	while read p1 p2 p3
	do
	  echo $p1$'\t'$p2$'\t'0$'\t'$p3$'\t'0$'\t'0$'\t'0 >> dU_velo_$number.txt 
	done < dU_re.txt
 	rm -f dU_re.txt dU.grd dU_re.grd dU_$number.txt tmp_dU.txt

	gmt psxy ciloto_prone_landslide.xy -JX10i -R$region-K -Wthinner,red -B0.05wesn --MAP_FRAME_TYPE=inside > $path/$psfile"_"$number.ps
	gmt psclip ciloto_prone_landslide.xy -J -R -O -K -V >> $path/$psfile"_"$number.ps
	gmt psvelo $path/dU_velo_$number.txt -R -J -W0.5p,red -Se0.2/0/6 -V -O -K >> $path/$psfile"_"$number.ps
	gmt psclip -J -R -C -O -V >> $path/$psfile"_"$number.ps
	gmt psconvert $psfile"_"$number.ps -Tg -W+k+t"vertical disp"+l256/-1 -E526 -V
	rm $path/$psfile"_"$number.ps
	mv dU_velo_$number.txt $path_ts/vertical/dU_velo/.
	mv $path/$psfile"_"$number.png $path_ts/vertical/.
	mv $path/$psfile"_"$number.kml $path_ts/vertical/.

	((i++)) 
   done

   convert -resize 50% -delay 20 -loop 0 $path_ts/vertical/dU_{001..020}.png dU.gif

else

	# create dU_velo.txt
	rm -f dU_velo.txt
	gmt xyz2grd dU.txt -Dm/m/m/1/0 -GdU.grd -I0.000046296 -R$region
	gmt grdsample dU.grd -I0.000138889 -GdU_re.grd
	gmt grd2xyz -R106.99/107.02/-6.73/-6.7 dU_re.grd > dU_re.txt
	while read p1 p2 p3
	do
	echo $p1$'\t'$p2$'\t'0$'\t'$p3$'\t'0$'\t'0$'\t'0 >> dU_velo.txt 
	done < dU_re.txt
	rm dU.grd dU_re.grd dU_re.txt

	# create map velocity for the vertical vector
	gmt psxy ciloto_prone_landslide.xy -JX10i -R$region -K -Wthinner,red -B0.05wesn --MAP_FRAME_TYPE=inside > $path/$psfile.ps
	gmt psclip ciloto_prone_landslide.xy -J -R -O -K -V >> $path/$psfile.ps
	gmt psvelo $path/dU_velo.txt -R -J -W0.5p,red -Se0.05/0/6 -V -O -K >> $path/$psfile.ps
	gmt psclip -J -R -C -O -V >> $path/$psfile.ps

	# convert ps to kml file
	cd $path
	gmt psconvert $psfile.ps -TG -W+k+t"vertical disp"+l256/-1 -E526 -V  #-Tg for white PNG, -TG for tranparent
	cd $dir

	##gmt psxy ciloto_prone_landslide.xy -JX10d -R -K -Wthinner,red -Bwesn > $psfile_gps
	##gmt psvelo gps_all_ver_geo.txt -R -J -W0.5p,black -Se0.25/0/6 -V -O >> $psfile_gps
	##gmt psconvert $psfile_gps -Tg -W+k+t"gps vertical disp"+l256/-1 -V

fi
