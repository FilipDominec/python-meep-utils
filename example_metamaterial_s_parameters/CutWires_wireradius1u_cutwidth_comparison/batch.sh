#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
model=SphereArray
cellsize=100e-6
thz=1e12
par='model=SphereArray resolution=2u simtime=100p radius=0u'

if [ -z "$skipsimulation" ]; then 
	for wc in 2 4 8 16 32 48; do
		mpirun -np $NP  ../../scatter.py $par wirethick=1u wirecut=${wc}e-6
		../../effparam.py
	done
fi

sharedoptions='effparam/*.dat --paramname wirecut --figsizey 2 --xeval x/1e12 --parameval param*1e6 --xlim1 0.4 --xlim2 1.2'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' --y2eval '0-y2' --ycol2 'imag N' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' --y2eval '0-y2' --ycol2 'imag eps' --ylim1 -2  --ylim2 5 \
	--paramlabel 'none' \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}$' --output ${PWD##*/}_eps.pdf
../../plot_multiline.py effparam/*wirecut=2.000e-06*dat --paramname wirecut --figsizey 2 --xeval x/1e12 --ylim1 0 --parameval param*1e6 \
	--xlabel "Frequency (THz)" --ycol 'real eps' --y2eval '0-y2' --ycol2 'imag eps' --ylim1 -2  --ylim2 5 \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}$' --output ${PWD##*/}_eps_r2.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' --y2eval '0-y2' --ycol2 'imag mu' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Permeability $\mu_{\text{eff}}$' --output ${PWD##*/}_mu.pdf

