#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
thz=1e12
cellsize=100e-6
par="resolution=2u radius=0 cellsize=$cellsize simtime=100p wirecut=2u"
#if [ -z "$skipsimulation" ]; then 
	#for wirethick in 1 2 4 8 16; do
		#mpirun -np $NP ../../scatter.py model=SphereWire $par wirethick=${wirethick}e-6 
		#../../effparam.py
	#done
#fi

sharedoptions="effparam/*.dat --paramname wirethick --parameval param*1e6 --figsizey 2 --xeval x/1e12"

../../plot_multiline.py $sharedoptions --paramlabel '$\rho_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf
../../plot_multiline.py $sharedoptions --paramlabel '$\rho_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf 


../../plot_multiline.py $sharedoptions --paramlabel '$\rho_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real N' --y2eval '0-y2' --ycol2 'imag N'  \
   	--ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz" 
../../plot_multiline.py $sharedoptions --paramlabel '$\rho_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real N' \
   	--ylabel 'Refractive index $N_{\text{eff}^{\prime}}$' --output ${PWD##*/}_nr.pdf \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz" 
../../plot_multiline.py $sharedoptions --paramlabel '$\rho_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}^{\prime\prime}}$' --output ${PWD##*/}_ni.pdf 


../../plot_multiline.py $sharedoptions --paramlabel '$\rho_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real eps' --y2eval '0-y2' --ycol2 'imag eps' --ylim1 -5 --ylim2 5  \
   	--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}$' --output ${PWD##*/}_eps.pdf
../../plot_multiline.py effparam/*thick=1.000e-06*.dat --paramname wirethick --parameval param*1e6 --figsizey 2 --xeval x/1e12 \
	--paramlabel '$\rho_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real eps' --y2eval '0-y2' --ycol2 'imag eps' --ylim1 -5 --ylim2 5  \
   	--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}$' --output ${PWD##*/}_eps_wt1.pdf
../../plot_multiline.py $sharedoptions --paramlabel '$\rho_w$ = %.1f $\upmu$m' --xlabel "Frequency (THz)" --ycol 'real mu' --y2eval '0-y2' --ycol2 'imag mu' --ylim1 -12 --ylim2 4  \
   	--ylabel 'Effective permeability $\mu_{\text{eff}}$' --output ${PWD##*/}_mu.pdf 
