#!/bin/bash
# 01.02.2018 NI create the script
# source: http://gmt.soest.hawaii.edu/boards/6/topics/3271 for local incidence and heading angles
# the local incidence angle will have greather value than the range of incidence angle (loc inc= ~50 deg)
# the local incidence angle is considered the topography and earth curvature variation (source: DEM)

# run script on 2D/3D folder
################ SET PARAMETER ###########################
input_asc=los_asc.txt
input_dsc=los_dsc.txt
ts_path=timeseries
aspect=aspect_val.txt
resolution=0.000046296 #0.000277778 #6.9445e-05
radius=0.000277778 #0.000046296 #0.000277778= 30 m radius
region=106.99/107.02/-6.73/-6.7 #-R106.97/107.02/-6.73/-6.7
#--------------------------------------------------------#
topo_asc=/home/isya/APPS/ciloto/Sentinel1/batch_asc/topo
topo_dsc=/home/isya/APPS/ciloto/Sentinel1/batch_dsc/topo
##########################################################

  dir=$(pwd)
  if [[ $# -ne 1 ]]; then
    echo ""
    echo "Usage: prepare_disp3d.sh #mode"
    echo ""
    echo "  script to prepare directory and process SBAS"
    echo ""
    echo "  example : prepare_disp3d.sh 1"
    echo "  or"
    echo "  example : prepare_disp3d.sh 2"
    echo ""
    echo "  Mode: 1 Prepare LOS asc and dsc files"
    echo "        2 Compute azimuth and incidence angle from the master scene asc & dsc"
    echo "        3 Convert and display the result on GMT environment"
    echo "        4 Prepare LOS asc and dsc files for time series (TS)"
    echo ""
    exit 1
  fi

mode=$1
echo "Mode -->" $mode

# going to Mode 1 
if [ $mode -eq 1 ]; then
echo "--------------> Prepare LOS asc and dsc files"
#gmt nearneighbor $input_asc -R$region -I$resolution -S$radius -Glos_asc.grd -N1
#gmt nearneighbor $input_dsc -R$region -I$resolution -S$radius -Glos_dsc.grd -N1
gmt surface $input_asc -R$region -I$resolution -Glos_asc_surface.grd -T0.25 -C0.1
gmt surface $input_dsc -R$region -I$resolution -Glos_dsc_surface.grd -T0.25 -C0.1
gmt grd2xyz -R$region los_asc_surface.grd > los_asc_surface.xyz
gmt grd2xyz -R$region los_dsc_surface.grd > los_dsc_surface.xyz
gmt xyz2grd $aspect -Ddegree/degree/degree/1/0 -Gaspect.grd -I0.001 -R$region -V
gmt grdsample aspect.grd -I$resolution -Gaspect_re.grd
gmt grd2xyz -R$region aspect_re.grd > aspect.xyz
#gmt grd2xyz -R$region los_asc.grd > los_asc_fin.xyz
#gmt grd2xyz -R$region los_dsc.grd > los_dsc_fin.xyz
# using azimuth from dem calculated by GMT
#gmt grdgradient dem.grd -Dno -Gdem_az.grd -R106.99/107.02/-6.73/-6.7 -V
#gmt grdsample dem_az.grd -I0.000046296 -Gdem_az_re.grd
#gmt grd2xyz -R$region dem_dir_re.grd > dem_az_fin.xyz
fi

# going to Mode 2 
if [ $mode -eq 2 ]; then
echo "--------------> Compute azimuth and incidence angle from the master scene"

# run on topo folder for ascending data
cd $topo_asc
gmt grdcut dem.grd -Gdem_cut.grd -R$region -V
gmt grdsample dem_cut.grd -I$resolution -Gdem_cut_fix.grd
gmt grd2xyz dem_cut_fix.grd > dem.xyz
SAT_look master.PRM < dem.xyz > dem_look.lltn
awk '{print $5}' dem_look.lltn > look_N
awk '{print $4}' dem_look.lltn > look_E
awk '{print $6}' dem_look.lltn > look_U

# calculate azimuth angle
gmt math look_N look_E ATAN2 = theta_look #(in radians)

# calculate incidence angle
gmt math look_E SQR = look_E_2
gmt math look_N SQR = look_N_2
gmt math look_E_2 look_N_2 ADD = look_E_N_2
gmt math look_E_N_2 SQRT = look_E_N_sqrt
gmt math look_U look_E_N_sqrt DIV = look_U_en
gmt math look_U_en ATAN = inc_angle #(in radians)
gmt math theta_look R2D = theta_look_degree
gmt math 180 theta_look_degree ADD -1 MUL = theta_look_degree_c
gmt math inc_angle R2D = inc_angle_degree

awk '{print $1,$2}' dem_look.lltn > ll
paste -d\  ll inc_angle_degree > inc_angle.lld
paste -d\  ll theta_look_degree_c > az_angle.lld
cp inc_angle.lld $dir/inc_angle_asc.lld
cp az_angle.lld $dir/az_angle_asc.lld

# inc and az angle convert to grid
gmt xyz2grd inc_angle.lld -Dm/m/m/1/0 -Ginc_angle_asc.grd -I$resolution -R$region -V
gmt xyz2grd az_angle.lld -Dm/m/m/1/0 -Gaz_angle_asc.grd -I$resolution -R$region -V

rm -f look_N look_E look_U look_E_2 look_N_2 look_E_N_2 look_E_N_sqrt look_U_en inc_angle inc_angle_degree theta_look_degree theta_look_degree_c ll dem.xyz dem_cut.grd theta_look dem_cut_fix.grd #dem_look.lltn

# run on topo folder for descending data
cd $topo_dsc
gmt grdcut dem.grd -Gdem_cut.grd -R$region -V
gmt grdsample dem_cut.grd -I$resolution -Gdem_cut_fix.grd
gmt grd2xyz dem_cut_fix.grd > dem.xyz
SAT_look master.PRM < dem.xyz > dem_look.lltn
awk '{print $5}' dem_look.lltn > look_N
awk '{print $4}' dem_look.lltn > look_E
awk '{print $6}' dem_look.lltn > look_U

# calculate azimuth angle
gmt math look_N look_E ATAN2 = theta_look #(in radians)

# calculate incidence angle
gmt math look_E SQR = look_E_2
gmt math look_N SQR = look_N_2
gmt math look_E_2 look_N_2 ADD = look_E_N_2
gmt math look_E_N_2 SQRT = look_E_N_sqrt
gmt math look_U look_E_N_sqrt DIV = look_U_en
gmt math look_U_en ATAN = inc_angle #(in radians)
gmt math theta_look R2D = theta_look_degree
gmt math 180 theta_look_degree ADD -1 MUL = theta_look_degree_c
gmt math inc_angle R2D = inc_angle_degree

awk '{print $1,$2}' dem_look.lltn > ll
paste -d\  ll inc_angle_degree > inc_angle.lld
paste -d\  ll theta_look_degree_c > az_angle.lld
cp inc_angle.lld $dir/inc_angle_dsc.lld
cp az_angle.lld $dir/az_angle_dsc.lld

# inc and az angle convert to grid
gmt xyz2grd inc_angle.lld -Dm/m/m/1/0 -Ginc_angle_dsc.grd -I$resolution -R$region -V
gmt xyz2grd az_angle.lld -Dm/m/m/1/0 -Gaz_angle_dsc.grd -I$resolution -R$region -V

rm -f look_N look_E look_U look_E_2 look_N_2 look_E_N_2 look_E_N_sqrt look_U_en inc_angle inc_angle_degree theta_look_degree theta_look_degree_c ll dem.xyz dem_cut.grd theta_look dem_cut_fix.grd #dem_look.lltn
 
cd $dir
fi

if [ $mode -eq 3 ]; then
echo "--------------> Convert and display the result on GMT environment"

gmt xyz2grd dU.txt -Dm/m/m/1/0 -GdU.grd -I$resolution -R$region -V
gmt xyz2grd dE.txt -Dm/m/m/1/0 -GdE.grd -I$resolution -R$region -V

#make cpt for dU
gmt grdinfo dU.grd > tmp
z_min=$(grep 'z_min:' tmp | awk '{print $3}')
z_max=$(grep 'z_max:' tmp | awk '{print $5}')
#gmt makecpt -Cpolar -T$z_min/$z_max/0.5 > dU.cpt
gmt makecpt -Cpolar -T-8/8/1 > dU.cpt
rm tmp
#make cpt for dE
gmt grdinfo dE.grd > tmp
z_min=$(grep 'z_min:' tmp | awk '{print $3}')
z_max=$(grep 'z_max:' tmp | awk '{print $5}')
#gmt makecpt -Cpolar -T$z_min/$z_max/0.5 > dE.cpt
gmt makecpt -Cpolar -T-5/5/1 > dE.cpt
rm tmp

#display dU and dE
gmt grdimage dU.grd -CdU.cpt -JX6i -I$resolution -B0.01WeSn -K -Y4 > dU.ps
gmt psscale -D2.3c/-1.2c/8c/0.4h -CdU.cpt -B4:"vertical displacament (mm)": -J -O -X3 -Y0 >> dU.ps
grd2kml.csh dU dU.cpt
gmt grdimage dE.grd -CdE.cpt -JX6i -I$resolution -B0.01WeSn -K -Y4 > dE.ps
gmt psscale -D2.3c/-1.2c/8c/0.4h -CdE.cpt -B2:"west-east displacament (mm)": -J -O -X3 -Y0 >> dE.ps
grd2kml.csh dE dE.cpt

fi

if [ $mode -eq 4 ]; then
echo "--------------> Prepare LOS asc and dsc files for time series (TS)"
cd $ts_path
#mkdir -p TXT
n=$(head -n 1 ts_asc.txt | awk '{ print NF}')
rm -f asc_surface.xyz dsc_surface.xyz
 
for (( i=1; i<=$n; i++ ))
do
    awk -v c=$i '{print $c}' ts_asc.txt > tmp_asc
    awk -v c=$i '{print $c}' ts_dsc.txt > tmp_dsc
    paste -d\,  lonlat_asc.txt tmp_asc > asc.in
    paste -d\,  lonlat_dsc.txt tmp_dsc > dsc.in
	gmt surface asc.in -R$region -I$resolution -Gasc_surface.grd -T0.1 -C0.5
	gmt surface dsc.in -R$region -I$resolution -Gdsc_surface.grd -T0.1 -C0.5
        if [ $i == 1 ]; then
	   gmt grd2xyz -R$region asc_surface.grd -Z > asc_surface.xyz
	   gmt grd2xyz -R$region dsc_surface.grd -Z > dsc_surface.xyz
        else
           gmt grd2xyz -R$region asc_surface.grd -Z > asc_surface_tmp.xyz
	   paste -d\,  asc_surface.xyz asc_surface_tmp.xyz > asc_surface_tmp2.xyz
           mv asc_surface_tmp2.xyz asc_surface.xyz
           gmt grd2xyz -R$region dsc_surface.grd -Z > dsc_surface_tmp.xyz
           paste -d\,  dsc_surface.xyz dsc_surface_tmp.xyz > dsc_surface_tmp2.xyz
           mv dsc_surface_tmp2.xyz dsc_surface.xyz 
        fi
    rm tmp_asc asc.in tmp_dsc dsc.in asc_surface_tmp.xyz dsc_surface_tmp.xyz #asc_surface.grd dsc_surface.grd
done 

cd ..
fi
