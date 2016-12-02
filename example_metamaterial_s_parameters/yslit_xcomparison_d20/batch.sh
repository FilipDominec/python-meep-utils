#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
thz=1e12
cellsize=100e-6
staticpar=(model=Fishnet resolution=2u simtime=350p cellsize=100u cellsizexy=100u slabcdist=0u slabthick=20e-6 ) 

if [ -z "$skipsimulation" ]; then 
	for	xhs in 4 10 20 40 80; do
		mpirun -np $NP ../../scatter.py "${staticpar[@]}" yholesize=inf xholesize=${xhs}u
		../../effparam.py --numstabcoef .97
	done
fi

sharedoptions=(effparam/*.dat --paramname xholesize --paramlabel '$d_x$ = %.1f $\upmu$m' --parameval param*1e6 --figsizey 2 
		--xlabel "Frequency (THz)" --xeval x/1e12 --xlim1 ' -0.5' --xlim2 3.0 --ylim1 0. --ylim2 1.)

../../plot_multiline.py "${sharedoptions[@]}"  --ycol '|r|' --ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf
../../plot_multiline.py "${sharedoptions[@]}"  --ycol '|t|' --ylabel 'Transmittance   $|t|$' --output ${PWD##*/}_t.pdf
