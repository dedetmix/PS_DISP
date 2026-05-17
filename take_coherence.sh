# 09.07.2018	NI	Do GMTSAR Interferogram TOPS Processing
#			Take Coherence files
#			Calculate the mean coherence

list=$1 
#input : sbas.list
raw=$2 
#raw path : /home/isya/APPS/ciloto/Sentinel1/batch_dsc/raw
mode=$3
# mode 1 to copy coherence from GMTSAR processing
# mode 2 to calculate the mean coherence from sbas.list
region=106.99/107.02/-6.73/-6.7

# example, type --> $ take_coherence.sh sbas.list /home/isya/APPS/ciloto/Sentinel1/batch_dsc/raw 1

if [ $mode -eq 1 ]; then
echo ""
echo "== taking coherence ... =="
mkdir -p corr

shopt -s extglob
IFS=" "
num=1
while read master slave
do

master_id=$(grep SC_clock_start $raw/$master.PRM | awk '{printf("%d",int($3))}')
slave_id=$(grep SC_clock_start $raw/$slave.PRM | awk '{printf("%d",int($3))}')

echo $master":"$slave > intf_exec.in  
intf_tops.csh intf_exec.in batch_tops.config
#cp intf_all/$master_id"_"$slave_id/corr.ps corr/corr_$master_id"_"$slave_id.ps
cp intf_all/$master_id"_"$slave_id/corr.grd corr/corr_$master_id"_"$slave_id.grd
cp intf_all/$master_id"_"$slave_id/corr_ll.grd corr/corr_$master_id"_"$slave_id"_ll".grd
#cp intf_all/$master_id"_"$slave_id/corr_ll.kml corr/corr_$master_id"_"$slave_id"_ll".kml
#cp intf_all/$master_id"_"$slave_id/corr_ll.png corr/corr_$master_id"_"$slave_id"_ll".png

rm -r intf_all/$master_id"_"$slave_id

done < $list

fi

if [ $mode -eq 2 ]; then
echo ""
echo "== calculate mean coherence ... =="

shopt -s extglob
IFS=" "
num=1
while read master slave
do

master_id=$(grep SC_clock_start $raw/$master.PRM | awk '{printf("%d",int($3))}')
slave_id=$(grep SC_clock_start $raw/$slave.PRM | awk '{printf("%d",int($3))}')
echo "$num Process $master_id and $slave_id correlation file"

  # crop correlation file
  gmt grdcut corr_$master_id"_"$slave_id"_ll".grd -R$region -Gcorr_cut.grd

  if [ $num == 1 ]; then
    gmt grdmath corr_cut.grd = sum.grd  
  else
    gmt grdmath sum.grd corr_cut.grd ADD = sumtmp.grd 
    mv sumtmp.grd sum.grd 
  fi
  (( num++ ))
  rm corr_cut.grd

done < $list

  (( num-- )) 
  gmt grdmath sum.grd $num DIV = mean_corr_ll.grd
  
  # plot the mean coherence on google earth
  gmt makecpt -T0./1/0.1 -Cgray -Z -N > corr.cpt
  echo "N  255   255   254" >> corr.cpt
  grd2kml.csh mean_corr_ll corr.cpt

fi
