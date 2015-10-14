#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
thz=1e12
mpirun -np $NP ../../scatter.py model=SphereArray resolution=5u radius=0 wirethick=4u
../../effparam.py
mpirun -np $NP ../../scatter.py model=SphereArray resolution=5u radius=0 wirethick=9u
../../effparam.py

sharedoptions='effparam/*.dat --paramname epsilon --paramlabel '\$\\varepsilon_r\$' --figsizey 2 --xeval x/1e12 --ylim1 0'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf # --color RdYlBu

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf  #--color RdYlBu_r

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf #--color PiYG_r


