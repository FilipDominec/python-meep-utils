#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
model=SphereArray
cellsize=100e-6
thz=1e12
par='model=SphereArray resolution=2u simtime=100p resolution=2u radius=0u'

for wc in 2 4 8 16 32 48; do
	mpirun -np $NP  ../../scatter.py $par wirethick=1u wirecut=${wc}e-6
	../../effparam.py
done

sharedoptions='effparam/*.dat --paramname wirecut --figsizey 2 --xeval x/1e12 --ylim1 0 --parameval param*1e6'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' --ycol2 'imag N' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$, $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_n.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' --ycol2 'imag eps' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime}$, $\varepsilon_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_eps.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' --ycol2 'imag mu' \
	--paramlabel '$d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Permeability $\mu_{\text{eff}}^{\prime}$, $\mu_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_mu.pdf

