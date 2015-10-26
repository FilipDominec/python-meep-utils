#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
model=SphereArray
cellsize=300e-6
thz=1e12
if [ -z "$skipsimulation" ]; then 
	for ff in `seq 0 .05 1 | tr , .`; do
		mpirun -np $NP ../../scatter.py model=Slab resolution=1u simtime=50p cellsize=$cellsize padding=10e-6 epsilon=4 fillfraction=$ff 
		../../effparam.py
	done
fi

sharedoptions="effparam/*.dat --paramname fillfraction --contours yes --colormap gist_earth_r --figsizex 4 --figsizex 4 --xeval x/1e12"

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf \
	--paramlabel 'Dielectric fill fraction $d_2/a$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance   $|t|$' --output ${PWD##*/}_t.pdf \
	--paramlabel 'Dielectric fill fraction $d_2/a$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"


# TODO
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf \
	--paramlabel 'Dielectric fill fraction $d_2/a$' \
	--overlayplot "(c/($cellsize*0.15*2*x*$thz))**2,(2*c/($cellsize*0.15*2*x*$thz))**2"
