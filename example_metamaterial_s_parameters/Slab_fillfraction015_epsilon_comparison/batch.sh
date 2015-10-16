#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
model=SphereArray
cellsize=300e-6
thz=1e12
for epsilon in 4 12 20; do
	mpirun -np $NP ../../scatter.py model=Slab fillfraction=.15 resolution=3u simtime=50p cellsize=$cellsize padding=100e-6 epsilon=4
	../../effparam.py
done

sharedoptions='effparam/*.dat --paramname epsilon --figsizey 2 --xeval x/1e12 --ylim1 0'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
	--paramlabel '$\varepsilon_r$' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '$\varepsilon_r$' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
	--paramlabel '$\varepsilon_r$' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag N' \
	--paramlabel '$\varepsilon_r$' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf

