########################################################

#Comparison of RV-TDT commands

########################################################

#1. Latino command that didn't work

#A. Latino - seg fault

"/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT" test1 -G latino.map -P latino.ped -M latino.map --adapt 500 --alpha 1e-05 --permut 2000 --lower_cutoff 0 --upper_cutoff 100 --minVariants 3 --maxMissRatio 1

#A. Euro - seg fault

"/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT" test1 -G euro.map -P euro.ped -M euro.map --adapt 500 --alpha 1e-05 --permut 2000 --lower_cutoff 0 --upper_cutoff 100 --minVariants 3 --maxMissRatio 1

########################################################

#2. Euro command that did work

#A. Euro - WORKS

#Includes './' in the input filename
'/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT' ./NA -G './euro.tped' -P './euro.ped' -M './euro.map' --adapt 500 --alpha 1e-05 --permut 2000 --lower_cutoff 0 --upper_cutoff 100 --minVariants 3 --maxMissRatio 1

#B. Latino - seg fault
'/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT' ./NA -G './latino.tped' -P './latino.ped' -M './latino.map' --adapt 500 --alpha 1e-05 --permut 2000 --lower_cutoff 0 --upper_cutoff 100 --minVariants 3 --maxMissRatio 1

########################################################

#3. Euro command that did work

#A. Euro - WORKS

#Gene name is test1
#Includes './' in the input filename
'/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT' ./test1 -G './euro.tped' -P './euro.ped' -M './euro.map' --adapt 500 --alpha 1e-05 --permut 2000 --lower_cutoff 0 --upper_cutoff 100 --minVariants 3 --maxMissRatio 1

########################################################