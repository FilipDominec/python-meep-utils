#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
thz=1e12
cellsize=100e-6
par="resolution=2u radius=0 cellsize=$cellsize simtime=20p"
for wirethick in 2 4 6 8; do
	mpirun -np $NP ../../scatter.py model=SphereWire $par wirethick=${wirethick}e-6 
	../../effparam.py
done

sharedoptions="effparam/*.dat --paramname wirethick --parameval 'param*1e6' --figsizey 2 --xeval x/1e12"

../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f' --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf # --color RdYlBu

../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f' --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf  #--color RdYlBu_r

../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f' --xlabel "Frequency (THz)" --ycol 'real N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
#
../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f' --xlabel "Frequency (THz)" --ycol 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf #--color PiYG_r


