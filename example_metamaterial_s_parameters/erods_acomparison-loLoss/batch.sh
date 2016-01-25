#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
thz=1e12
cellsize=100e-6
par="resolution=1u simtime=150p loss=.1"
if [ -z "$skipsimulation" ]; then 
	for a in 80 90 120; do
		mpirun -np $NP ../../scatter.py $par model=RodArray orientation=E radius=10e-6  cellsize=${a}u 
		../../effparam.py
	done
fi

sharedoptions="effparam/*.dat --paramname cellsize  --parameval param*1e6 --figsizey 2 --xeval x/1e12"

../../plot_multiline.py $sharedoptions --paramlabel '$a$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf
../../plot_multiline.py $sharedoptions --paramlabel '$a$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf 


../../plot_multiline.py $sharedoptions --paramlabel '$a$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real N' --y2eval '0-y2' --ycol2 'imag N' \
--ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf 
    #--overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz" 

../../plot_multiline.py $sharedoptions --paramlabel '$a$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real eps' --y2eval '0-y2' --ycol2 'imag eps' --ylim1 -10 --ylim2 10  \
   	--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_eps.pdf \
    --overlayplot "1-.68**2/x**2"  
../../plot_multiline.py $sharedoptions --paramlabel '$a$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real eps' --ylim1 -12 --ylim2 3  \
   	--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_epsr.pdf \
    --overlayplot "1-.68**2/x**2"  
../../plot_multiline.py $sharedoptions --paramlabel '$a$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag eps' --ylim1 -12 --ylim2 3  \
   	--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_epsi.pdf \


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' --y2eval '0-y2' --ycol2 'imag mu' \
	--ylim1 -10 --ylim2 10 \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Permeability $\mu_{\text{eff}}$' --output ${PWD##*/}_mu.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' \
	--ylim1 -5 --ylim2 5 \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Permeability $\mu_{\text{eff}}^{\prime}$' --output ${PWD##*/}_mur.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag mu' \
	--ylim1 -5 --ylim2 5 \
	--paramlabel '$a$ = %.1f $\upmu$m' --ylabel 'Permeability $\mu_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_mui.pdf
