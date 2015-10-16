#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
cellsize=300e-6
thz=1e12
for eps in `seq 1 .5 25 | tr , .`; do
	mpirun -np $NP ../../scatter.py model=Slab resolution=1u simtime=50p cellsize=$cellsize padding=10e-6 fillfraction=.15 epsilon=$eps
	../../effparam.py
done

sharedoptions="effparam/*.dat --paramname epsilon --contours yes --colormap gist_earth_r --figsizex 4 --figsizex 4 --xeval x/1e12 --ylim1 0"

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf \
	--paramlabel 'Dielectric permittivity $\varepsilon_r$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance   $|t|$' --output ${PWD##*/}_t.pdf \
	--paramlabel 'Dielectric permittivity $\varepsilon_r$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf \
	--paramlabel 'Dielectric permittivity $\varepsilon_r$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"
