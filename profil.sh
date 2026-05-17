#!/bin/bash

### SET PARAMETER ###
region=106.9932/106.9947/-6.7075/-6.7055
path_ts=timeseries
ps=profil_puncak_pass.ps
number=087
#####################

# create tracking A-B
#gmt project -C106.993923111000000/-6.705602979010000 -E106.994241840000000/-6.707167117700000 -G0.001 -Q -N > track_puncakpass.txt
# track dU_surface and DEM
#awk {'print $1,$2'} track_puncakpass.txt | gmt grdtrack -Gciloto_30m_subset.grd -h > table_ver_elev.txt
#awk {'print $1,$2'} track_puncakpass.txt | gmt grdtrack -GdU_surface.grd -h > table_dU.txt

# edit at open office for table_ver_elev.txt table_ver_selisih.txt

gmt psxy area_puncakpass.txt -JM20 -B0.001 -P -Y3 -R$region -Wthinner,red -K > $ps
gmt psclip area_puncakpass.txt -J -B0.001 -P -R$region -O -K >> $ps
gmt grdimage $path_ts/velocity/GRD/dU_surface_$number.grd -J -Cver.cpt -B -P -R$region -O -K  >> $ps
gmt psvelo $path_ts/velocity/TXT/dEN_velo_$number.txt -R -J -W0.5p,black -Se0.05/0/6 -O -K >> $ps
gmt psclip -C -O -K >> $ps
gmt psxy -J -R -B -W0.4p,purple -O -K profil_puncakpass.txt >> $ps
gmt psxy -J -R -B -Sc0.01i -Gblack -O -K profil_puncakpass.txt >> $ps
gmt psscale -DjCB+w4i/1c+o0/2c+h -Cver.cpt -R -J -Bx5+l"cumulative vertical disp (mm)" -O -K >> $ps

awk {'print $4,$3'} table_ver_elev.txt | gmt psxy -R0/179/1420/1490 -Bx+l"Distance from ridge 170 (m)" -BWS -Byaf100+l"elevation (scale to mm)" \
-JX8.5i/4i -O -K -Wgray -Sb1.5c -Y11.5i >> $ps
awk {'print $1,$2'} table_ver_selisih.txt | gmt psxy -R -B \
-J -O -K -Wred >> $ps
awk {'print $4,$3'} table_ver_elev.txt | gmt psxy -R -B \
-J -O -Wpurple >> $ps
