#!/bin/bash

ps=3d_arrow
region=90.544296/90.557503/30.4398/30.4592
path=/home/isya/APPS/China/vonEike/timeseries/plot_3d_SF_opt/DATA

method=SF_opt
source_dE=/home/isya/APPS/China/vonEike/timeseries/dE_ts_3d_"$method".txt
source_dN=/home/isya/APPS/China/vonEike/timeseries/dN_ts_3d_"$method".txt
source_lonlat=/home/isya/APPS/China/vonEike/timeseries/lonlat.txt
low_elev=4500
top_elev=6000

ln -s $source_dE .
ln -s $source_dN .
ln -s $source_lonlat .

# create DEM
#gmt grdsample dem.grd -Gdem_re.grd -I0.000277778 -R90.544296/90.557503/30.4398/30.4592
#gmt grdgradient dem_re.grd -A1 -Gdem_intens.grd -Nt0.9 -fg

mkdir JPG PS

for i in {1..95}
do
   number=$(printf '%03d\n' $i)

   # plot 3D
   gmt grdview dem_re.grd -Idem_intens.grd -G$path/dU_surface_$number.grd -R$region/$low_elev/$top_elev -JM10 -JZ5i -Qc300 -B0.005 -BwSEnZ -Bz1000+l"Elevation (m)" -p165/30 -C$path/ver.cpt -N$low_elev+glightgray -K -P > "$ps"_$number.ps
   # plot velocity
   ./velo_legend.sh
   gmt psxyz velo_hor_$i.txt -R$region/$low_elev/$top_elev -JM -JZ -P -Sv0.1+e -p -O >> "$ps"_$number.ps
   # plot point on elevation
   #gmt grdtrack -Gdem_re.grd area_puncakpass.txt | gmt psxyz -JM -JZ -p -R$region/5000/6500 -W2p,red -O >> "$ps"_$number.ps

   psconvert -Tj -E256 "$ps"_$number.ps
   mv *.ps PS/.
   mv *.jpg JPG/.
   rm -f theta_$i.txt dR_$i.txt
done

rm -f velo_hor* dR* theta*

#make gif
#convert to gif using imagemagick
#	convert -delay 20 -loop 0 *.jpg tropo_delay_aguablanca.gif
# if it's killed due to memory, split to every year or month
#	$ mkdir JPG
#	$ mv .jpg JPG/.
#	$ mkdir JPG1
#	$ for i in {1..365}; do number=$(printf "%03d" "$i"); ln -s $PWD/JPG/era.txt_$i.jpg JPG1/era_$number.jpg; done
# or with certain intervals and crop
#       $ for i in {1..365..14}; do number=$(printf "%03d" "$i"); ln -s $PWD/JPG/era.txt_$i.jpg JPG1/era_$number.jpg; convert JPG1/era_$number.jpg -crop 1170x827-350+200 JPG1/era_crop_$number.jpg; done
#	$ convert -delay 20 -loop 0 JPG1/era_crop_*.jpg tropo_delay_aguablanca_01.gif
# do every year or month

# short the duration of GIF and CROP
# S gifsicle --info mexico.gif
# $ gifsicle -U mexico.gif `seq -f "#%g" 0 14 364` -O2 -o output.gif
# $ convert phase_delay_mexico.gif -coalesce -repage 0x0 -crop 1170x827-225+100 +repage output.gif
