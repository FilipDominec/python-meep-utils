#!/bin/bash
if [ -z $NP ] ; then NP=1 ; fi			 # number of processors
par='model=Slab  radius=20n resolution=2n cellsize=100n' ## TODO

for x in `seq TODO` ; dob
	mpirun -np $NP   python ../../scatter.py $par ## TODO
	../../effparam.py 
done

../../plot_multiline.py PlasmonicSpheres_*dat  \
        --ycol '|t|'  --xlabel 'Frequency [Hz]' --ylabel 'transmission amplitude $|s_{12}|$'  \
        --paramname cellsize  --paramlabel '$\Delta x$, $\Delta y$ = %.0f nm' --paramunit 1e-9
