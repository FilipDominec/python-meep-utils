#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
model=SphereArray
cellsize=50e-6
#mpirun -np $NP ../../scatter.py model=SphereArray resolution=4u simtime=50p wirethick=10u cellsize=$cellsize padding=0e-6 radius=13e-6
mpirun -np $NP ../../scatter.py model=Slab fillfraction=1 resolution=5u simtime=50p cellsize=$cellsize padding=100e-6 epsilon=49

../../effparam.py

sharedoptions='effparam/*.dat --paramname radius --paramlabel none --figsizey 2 --xeval x/1e12'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf --color RdYlBu

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf --color RdYlBu_r

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' --ylim2 5. \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf --color PiYG_r \
    --overlayplot "2.998e8/4./${cellsize}/(x*1e12)"

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag N' --ylim2 5. \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf --color PiYG_r


