#!/bin/bash

###################### Loop to take xyz, remove first two lines and make into input file #######

for files in *.xyz
do

        newfile=$(echo $files | sed 's/\.[^.]*$//')
        cp $files $newfile
echo '$end

$rem
EXCHANGE M062X
BASIS 6-31+g*
DFT_D EMPIRICAL_GRIMME3
JOBTYPE SP
$end' >> $newfile
        sed '1,2d' $newfile > $newfile.inp
        rm $newfile
done

##########################Generates input parameters, adds to first line #######################

sed -i '1s/^/-1 2\n/' *.inp
sed -i '1s/^/$molecule\n/' *.inp

################################ Adds to the end ###############################################
############################ Optimisation settings #############################################

