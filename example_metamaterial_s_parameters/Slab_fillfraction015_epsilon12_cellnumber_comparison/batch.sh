#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
model=SphereArray
cellsize=300e-6
thz=1e12
if [ -z "$skipsimulation" ]; then 
	for cellnumber in 16 8 4 2 1; do
		mpirun -np $NP ../../scatter.py model=Slab fillfraction=.15 resolution=3u simtime=250p cellsize=$cellsize cellnumber=$cellnumber padding=100e-6 epsilon=4
		../../effparam.py
	done
fi

sharedoptions='effparam/*.dat --paramname cellnumber --figsizey 2 --xeval x/1e12'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
	--paramlabel 'cell number $m = %d$' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf 
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel 'cell number $m = %d$' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf  


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
	--paramlabel 'cell number $m = %d$' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' \
	--paramlabel 'cell number $m = %d$' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf 

