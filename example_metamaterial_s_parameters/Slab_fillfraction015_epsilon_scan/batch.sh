#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
cellsize=300e-6
thz=1e12
if [ -z "$skipsimulation" ]; then 
	for eps in `seq 1 .5 25 | tr , .`; do
		mpirun -np $NP ../../scatter.py model=Slab resolution=1u simtime=50p cellsize=$cellsize padding=10e-6 fillfraction=.15 epsilon=$eps
		../../effparam.py
	done
fi

sharedoptions="effparam/*.dat --paramname epsilon --contours yes --colormap gist_earth_r --figsizex 4 --figsizex 4 --xeval x/1e12"

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf \
	--paramlabel 'Dielectric permittivity $\varepsilon_r$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance   $|t|$' --output ${PWD##*/}_t.pdf \
	--paramlabel 'Dielectric permittivity $\varepsilon_r$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' --y2eval '0-y' --ycol2 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf \
	--paramlabel 'Dielectric permittivity $\varepsilon_r$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime}$' --output ${PWD##*/}_ni.pdf \
	--paramlabel 'Dielectric permittivity $\varepsilon_r$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf \
	--paramlabel 'Dielectric permittivity $\varepsilon_r$'
