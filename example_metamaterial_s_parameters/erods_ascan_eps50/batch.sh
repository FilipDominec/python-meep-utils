#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
cellsize=100e-6
thz=1e12
par="resolution=2u simtime=200p loss=.1 epsilon=50"

if [ -z "$skipsimulation" ]; then 
	for a in `seq 20 4 200`; do
		mpirun -np $NP ../../scatter.py $par model=RodArray orientation=E radius=10e-6  cellsize=${a}u 
		../../effparam.py
	done
fi

sharedoptions='effparam/*.dat --paramname cellsize --contours yes --numcontours 25 --colormap gist_earth --figsizex 4 --figsizey 3 --parameval param*1e6 --xeval x/1e12'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf \
	--paramlabel 'Unit cell size $a$ ($\upmu$m)' --overlayplot "80,95,120"
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance   $|t|$' --output ${PWD##*/}_t.pdf \
	--paramlabel 'Unit cell size $a$ ($\upmu$m)' --overlayplot "80,95,120"


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime}$' --output ${PWD##*/}_nr.pdf \
	--paramlabel 'Unit cell size $a$ ($\upmu$m)' --overlayplot "80,95,120"

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf \
	--ylim1 " -1.5"  --ylim2 0.4 \
	--paramlabel 'Unit cell size $a$ ($\upmu$m)' --overlayplot "80,95,120"

