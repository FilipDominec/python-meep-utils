#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
cellsize=100e-6
thz=1e12
par="resolution=2u simtime=200p loss=1"

if [ -z "$skipsimulation" ]; then 
	for a in `seq 20 8 200`; do
		mpirun -np $NP ../../scatter.py $par model=RodArray orientation=E radius=10e-6  cellsize=${a}u 
		../../effparam.py
	done
fi

sharedoptions='effparam/*.dat --paramname cellsize --parameval param*1e6 --figsizey 2 --xeval x/1e12 --xlim1 .25 --xlim2 1.25'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' --y2eval '0-y2' --ycol2 'imag N' --ylim1 -5 --ylim2 5\
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' --ylim2 3 \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' --y2eval '0-y2' --ycol2 'imag eps' --ylim1 -5 --ylim2 5 \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Permittivity $\varepsilon_{\text{eff}}$' --output ${PWD##*/}_eps.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_epsr.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag eps' \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_epsi.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' --y2eval '0-y2' --ycol2 'imag mu' \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Permeability $\mu_{\text{eff}}$' --output ${PWD##*/}_mu.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Permeability $\mu_{\text{eff}}^{\prime}$' --output ${PWD##*/}_mur.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag mu' \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Permeability $\mu_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_mui.pdf
