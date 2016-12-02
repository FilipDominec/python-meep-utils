#!/bin/bash
if [ -z $NP ] ; then NP=1 ; fi             # number of processors
thz=1e12
cellsize=100e-6
staticpar=(model=Slab resolution=5u ) 

if [ -z "$skipsimulation" ]; then 
	mpirun -np $NP ../ldos.py "${staticpar[@]}" yholesize=inf xholesize=${xhs}u
fi

	#sharedoptions=(effparam/*.dat --paramname xholesize --paramlabel '$d_x$ = %.1f $\upmu$m' --parameval param*1e6 --figsizey 2 
		#--xlabel "Frequency (THz)" --xeval x/1e12 --xlim1 ' -0.5' --xlim2 3.0 --ylim1 0. --ylim2 1.)
#
#../../plot_multiline.py "${sharedoptions[@]}"  --ycol '|r|' --ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf
#../../plot_multiline.py "${sharedoptions[@]}"  --ycol '|t|' --ylabel 'Transmittance   $|t|$' --output ${PWD##*/}_t.pdf
