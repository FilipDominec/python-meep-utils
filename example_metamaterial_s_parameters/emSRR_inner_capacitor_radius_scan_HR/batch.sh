#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
cellsize=100e-6
thz=1e12
par="resolution=4u simtime=100p cellsize=$cellsize"

if [ -z "$skipsimulation" ]; then 
	for icr in `seq 50 10 200` `seq 55 10 200`; do
		mpirun -np $NP  ../../scatter.py $par  \
				model=ESRRArray comment="emcSRR" cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
				insplitting=6e-6 incapacitorr=${icr}e-7 wirethick=0 radius=40e-6 srrthick=10e-6
		../../effparam.py
	done
fi


sharedoptions='effparam/*.dat  --paramname incapacitorr --contours yes --numcontours 15 --colormap gist_earth --figsizex 4 --figsizey 3 --parameval param*1e6 --xeval x/1e12 --xlim2 1.5'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' --y2eval '0-y2' --ycol2 'imag N' --ylim1 -5 --ylim2 5\
	--paramlabel 'Inner capacitor radius $\rho_c$' --ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
	--paramlabel 'Inner capacitor radius $\rho_c$' --ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' --ylim1 -2 --ylim2 -0.001 \
	--paramlabel 'Inner capacitor radius $\rho_c$' --ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf

