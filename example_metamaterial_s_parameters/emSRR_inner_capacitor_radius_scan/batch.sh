#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
cellsize=100e-6
thz=1e12
par="resolution=4u simtime=200p cellsize=$cellsize"

if [ -z "$skipsimulation" ]; then 
	for icr in 6 8 10 18; do
		mpirun -np $NP  ../../scatter.py $par  \
				model=ESRRArray comment="emcSRR" cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
				insplitting=6e-6 incapacitorr=${icr}e-6 wirethick=0 radius=40e-6 srrthick=10e-6
		../../effparam.py
	done
fi


sharedoptions='effparam/*.dat --paramname incapacitorr --parameval param*1e6 --figsizey 2 --xeval x/1e12 --xlim1 .25 --xlim2 1.25'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' --y2eval '0-y2' --ycol2 'imag N' --ylim1 -5 --ylim2 5\
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' --ylim2 3 \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' --y2eval '0-y2' --ycol2 'imag eps' --ylim1 -5 --ylim2 5 \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Permittivity $\varepsilon_{\text{eff}}$' --output ${PWD##*/}_eps.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_epsr.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag eps' \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_epsi.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' --y2eval '0-y2' --ycol2 'imag mu' \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Permeability $\mu_{\text{eff}}$' --output ${PWD##*/}_mu.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Permeability $\mu_{\text{eff}}^{\prime}$' --output ${PWD##*/}_mur.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag mu' \
	--paramlabel '$\rho_c$ = %.1f $\upmu$m' --ylabel 'Permeability $\mu_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_mui.pdf
