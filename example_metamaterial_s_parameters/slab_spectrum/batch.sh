#!/bin/bash
if [ -z $NP ] ; then NP=1 ; fi
#mpirun -np $NP  ../../scatter.py model=Slab cellsize=500u
../../effparam.py

../../plot_multiline.py S*dat  \
        --ycol '|t|'  --xlabel 'Frequency [Hz]' --ylabel 'transmission amplitude $|s_{12}|$'  \
        --paramname cellsize  --paramlabel '$\Delta x$, $\Delta y$ = %.0f nm' --paramunit 1e-9
