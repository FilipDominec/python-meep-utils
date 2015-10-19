#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
thz=1e12
cellsize=100e-6
par="resolution=1u radius=0 cellsize=$cellsize simtime=50p"
for wirethick in 1 2 4 8 16; do
	mpirun -np $NP ../../scatter.py model=SphereWire $par wirethick=${wirethick}e-6 
	../../effparam.py
done

sharedoptions="effparam/*.dat --paramname wirethick --parameval param*1e6 --figsizey 2 --xeval x/1e12"

../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf

../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf 

../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz" --ylim2 1.5  
#
../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'imag N' --ylim1 -0.5 --ylim2 2.5 \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf

../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real eps' --ylim1 -12 --ylim2 3  \
   	--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_epsr.pdf \
    --overlayplot "1-.6**2/x**2"  

../../plot_multiline.py $sharedoptions --paramlabel '$r_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'imag eps' --ylim1 -12 --ylim2 3  \
   	--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_epsi.pdf \

