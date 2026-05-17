#!/bin/bash

# Modify the number values for dE and dN plot in GMT and Google Earth (time series data)

number=$1
region=$2
resolution=0.000138889

## dE
# resample
mv dE_$number.txt dE_"$number"_ORI.txt
gmt xyz2grd dE_"$number"_ORI.txt -Ddegree/degree/degree/1/0 -GdE_resamp.grd -I$resolution -R$region -V
gmt grd2xyz -R$region dE_resamp.grd > dE_$number.txt
sed -i $'s/\t/,/g' dE_$number.txt
rm dE_resamp.grd dE_"$number"_ORI.txt
# filter big value
#rm -f dE_"$number"_edit.txt
#while read long lat value
#do
#   if [ $value == "NaN" ]; then
#      echo $long $lat $value >> dE_"$number"_edit.txt
#   else
#      if (( $(echo "$value < -300" |bc -l) || $(echo "$value > 300" |bc -l) ))
#      then
#         echo $long $lat "NaN" >> dE_"$number"_edit.txt
#      else 
#         echo $long $lat $value >> dE_"$number"_edit.txt
#      fi
#   fi
#done < dE_$number.txt
#mv dE_$number.txt dE_"$number"_bf_filter.txt
#mv dE_"$number"_edit.txt dE_$number.txt

## dN
# resample
mv dN_$number.txt dN_"$number"_ORI.txt
gmt xyz2grd dN_"$number"_ORI.txt -Ddegree/degree/degree/1/0 -GdN_resamp.grd -I$resolution -R$region -V
gmt grd2xyz -R$region dN_resamp.grd > dN_$number.txt
sed -i $'s/\t/,/g' dN_$number.txt
rm dN_resamp.grd dN_"$number"_ORI.txt
# filter big value
#rm -f dN_"$number"_edit.txt
#while read long lat value
#do
#   if [ $value == "NaN" ]; then
#      echo $long $lat $value >> dN_"$number"_edit.txt
#   else
#      if (( $(echo "$value < -300" |bc -l) || $(echo "$value > 300" |bc -l) ))
#      then
#         echo $long $lat "NaN" >> dN_"$number"_edit.txt
#      else 
#         echo $long $lat $value >> dN_"$number"_edit.txt
#      fi
#   fi
#done < dN_$number.txt
#mv dN_$number.txt dN_"$number"_bf_filter.txt
#mv dN_"$number"_edit.txt dN_$number.txt
