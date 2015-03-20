#!/bin/bash
# Note: in multiprocessing environment, run this script e.g. as `mpirun -np 2 ./batch.sh'

COMMAND='python ../scatter.py resolution=4u simtime=100p wirethick=10u'

## Generate frequency-domain results
for ff in `seq 1000 10 1300`; do  
     $COMMAND frequency=${ff}e9
done

## Extract the reflection (|s11|, second column) and transmission (|s12|, fourth column)
cat *frequency=*dat > all.dat
cat all.dat  |  sed -e '/#/d'  |  cut -d' ' -f'1,2'  |  sort -g  >  r.dat
cat all.dat  |  sed -e '/#/d'  |  cut -d' ' -f'1,4'  |  sort -g  >  t.dat

## Gather the frequency-domain E_x shapes
mkdir png
mv    SphereArray*frequency*/*png png/

## Clean up
rm    SphereArray*frequency*.dat
rm    SphereArray*frequency*.png
rm -r SphereArray*frequency*/

## Run one time-domain simulation for comparison
$COMMAND 
python ../effparam.py
