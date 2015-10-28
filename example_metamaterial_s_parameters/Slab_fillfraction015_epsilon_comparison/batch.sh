#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
model=SphereArray
cellsize=300e-6
thz=1e12
if [ -z "$skipsimulation" ]; then 
	for epsilon in 4 12 20; do
		mpirun -np $NP ../../scatter.py model=Slab fillfraction=.15 resolution=3u simtime=50p cellsize=$cellsize padding=100e-6 epsilon=$epsilon
		../../effparam.py
	done
fi

sharedoptions='effparam/*.dat --paramname epsilon --figsizey 2 --xeval x/1e12'

../../plot_multiline.py $sharedoptions --ycol '|r|' --xlabel "Frequency (THz)" \
	--paramlabel '$\varepsilon_r$' \
   	--ylabel 'Reflectance $|r|$' --output ${PWD##*/}_r.pdf
../../plot_multiline.py $sharedoptions --ycol '|t|' --xlabel "Frequency (THz)" \
	--paramlabel '$\varepsilon_r$' \
   	--ylabel 'Transmittance $|t|$' --output ${PWD##*/}_t.pdf


../../plot_multiline.py $sharedoptions --ycol 'real N' --ycol2 'imag N' --y2eval '0-y2' --xlabel "Frequency (THz)" \
	--paramlabel '$\varepsilon_r$' \
   	--ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  --figsizey 3 

