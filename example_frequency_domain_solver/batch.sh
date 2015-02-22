#!/bin/bash

COMMAND='mpirun -np 1   python ../scatter.py resolution=4u simtime=50p padding=00u'

## Generate frequency-domain results
for ff in `seq 1250 10 1300`; do  
     $COMMAND frequency=${ff}e9
done

## Extract the reflection (|s11|, second column) and transmission (|s12|, fourth column)
cat *frequency=*dat > all.dat
cat all.dat  |  sed -e '/#/d'  |  cut -d' ' -f'1,2'  |  sort -g  >  r.dat
cat all.dat  |  sed -e '/#/d'  |  cut -d' ' -f'1,4'  |  sort -g  >  t.dat

## Gather the frequency-domain E_x shapes
mkdir png
mv SphereArray_/*png png/

## Clean up
rm -r SphereArray*frequency*.dat

## Run one time-domain simulation for comparison
$COMMAND 
