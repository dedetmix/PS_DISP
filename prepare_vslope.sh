# GMT process to create clope for vslope computation

resolution=0.000046296
region=106.996974/107.011761/-6.720616/-6.710175

cd topo # either from batch_asc or batch_dsc
gmt grdgradient dem.grd -fg -Sslope_tmp.grd -D
gmt grdmath slope_asc_tmp.grd ATAN PI DIV 180 MUL = slope.grd
rm slope_tmp.grd

gmt grdsample slope.grd -I$resolution -R$region -Gslope_re.grd
gmt grd2xyz -R$region slope_re.grd > slope.xyz

# go to matlab
cd timeseries
# >> prepare_vce.m
 

