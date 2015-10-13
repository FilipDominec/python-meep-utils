#!/bin/bash
if [ -z $NP ] ; then NP=1 ; fi
for x in 1 2 3 4 5 10; do  
	mpirun -np $NP  ../../scatter.py model=Slab cellsize=500u cellnumber=$x
done
../../effparam.py

../../plot_multiline.py S*dat  \
        --ycol '|t|'  --xlabel 'Frequency [THz]' --xeval x/1e12 --ylabel 'transmission amplitude $|s_{12}|$'  \
        --paramname cellnumber  --paramlabel 'n_{cell} = %d' --paramunit 1
