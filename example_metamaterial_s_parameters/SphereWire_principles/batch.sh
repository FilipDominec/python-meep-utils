#!/bin/bash
if [ -z $NP ] ; then NP=1 ; fi			 # number of processors
par='model=SphereArray resolution=4u simtime=50p'

mpirun -np $NP  ../../scatter.py $par wirethick=4u radius=0u  comment=OnlyWire                   ;  ../../effparam.py
../effparam
mpirun -np $NP  ../../scatter.py $par wirethick=0u radius=13u comment=OnlySphere                 ;  ../../effparam.py
../effparam
mpirun -np $NP  ../../scatter.py $par wirethick=4u radius=13u comment=SphereWireNIM              ;  ../../effparam.py
../effparam


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf # --color RdYlBu

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf  #--color RdYlBu_r

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "-c/2/$cellsize/x/$thz, c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf #--color PiYG_r

