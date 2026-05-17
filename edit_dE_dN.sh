#!/bin/bash

# Modify the values for dE and dN plot in GMT and Google Earth

# dE
rm -f dE_3d_SF_edit.txt
while read long lat value
do
   if [ $value == "NaN" ]; then
      echo $long $lat $value >> dE_3d_SF_edit.txt
   else
      if (( $(echo "$value < -300" |bc -l) || $(echo "$value > 300" |bc -l) ))
      then
         echo $long $lat "NaN" >> dE_3d_SF_edit.txt
      else 
         echo $long $lat $value >> dE_3d_SF_edit.txt
      fi
   fi
done < dE_3d_SF.txt
mv dE_3d_SF.txt dE_3d_SF_bf_filter.txt
mv dE_3d_SF_edit.txt dE_3d_SF.txt

# dN
rm -f dN_3d_SF_edit.txt
while read long lat value
do
   if [ $value == "NaN" ]; then
      echo $long $lat $value >> dN_3d_SF_edit.txt
   else
      if (( $(echo "$value < -300" |bc -l) || $(echo "$value > 300" |bc -l) ))
      then
         echo $long $lat "NaN" >> dN_3d_SF_edit.txt
      else
         echo $long $lat $value >> dN_3d_SF_edit.txt
      fi
   fi
done < dN_3d_SF.txt
mv dN_3d_SF.txt dN_3d_SF_bf_filter.txt
mv dN_3d_SF_edit.txt dN_3d_SF.txt
