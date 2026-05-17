ps=3d_arrow.ps
region=106.9936/106.9948/-6.7075/-6.7055
number=095
i=95

# create DEM
#gmt grdsample dem.grd -Gdem_re.grd -I0.000046296 -R106.9933/106.9948/-6.7075/-6.7055
#gmt grdgradient dem_re.grd -A1 -Gdem_intens.grd -Nt0.9 -fg

# plot 3D
gmt grdview dem_re.grd -Idem_intens.grd -GdU_surface_$number.grd -R$region/1400/1500 -JM10 -JZ3.5i -Qc300 -B0.001 -BwsenZ1+n -Bz20+l"Elevation (m)" -p160/30 -Cver.cpt -N1400+glightgray -K -P > $ps
#gmt grdview dem_re.grd -Idem_intens.grd -R$region/1400/1500 -JM14 -JZ2i -Qsm -Bz20+l"elevation (m)" -p170/30 -CzeroShift.cpt -K -P -N1400+glightgray -B0.001 -BnESwZ > $ps

# plot velocity
#gmt psxyz velo_hor.txt -R$region/1400/1500 -B0.1 -BWeSnZ -Bz100 -JX10 -JZ10 -P -Sv0.25+e -p135/30 > test.ps
./velo_legend.sh
gmt psxyz velo_hor_$i.txt -R$region/1400/1500 -B -JM -JZ -P -Sv0.1+e -p -K -O >> $ps

# plot point on elevation
gmt grdtrack -Gdem_re.grd area_puncakpass.txt | gmt psxyz -JM -JZ -p -R$region/1400/1500 -W2p,red -O >> $ps

# plot 2D
#gmt psxy area_puncakpass.txt -JM8 -B0.001 -BwESn -p170/30 -R$region -Wthinner,red -O -K -X12 >> $ps
#gmt psclip area_puncakpass.txt -J -p -R$region -O -K >> $ps
#gmt grdimage dU_surface_$number.grd -J -Cver.cpt -p -R$region -O -K  >> $ps
#gmt psvelo dEN_velo_$number.txt -R -J -W0.05p,black -p -Se0.01/0/6 -O -K >> $ps
#gmt psclip -C -O -K >> $ps
#gmt psxy -J -R -W0.4p,purple -p -O -K profil_puncakpass.txt >> $ps
#gmt psxy -J -R -Sc0.01i -Gblack -p -O -K profil_puncakpass.txt >> $ps
#gmt psscale -DjCB+w4i/1c+o0/2c+h -Cver.cpt -R -J -Bx25+l"cumulative vertical disp (mm)" -p -O -K >> $ps
#echo 106.9939 -6.7070 50 0 0 0 0 50 mm | gmt psvelo -J -R -Se0.01/0/12 -W1p,black -L -V -p -K -O >> $ps

# plot legend
#echo "106.9937 -6.7071 Profil A-B" | gmt pstext -R -J -F+f12p,Helvetica,black+jLM -O -K >> $ps
#echo "106.99392 -6.7071 ____" | gmt pstext -R -J -F+f12p,Helvetica-Bold,purple+jLM -O >> $ps

psconvert -Tj -E256 $ps
