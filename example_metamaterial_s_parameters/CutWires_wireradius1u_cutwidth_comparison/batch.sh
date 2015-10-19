#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
model=SphereArray
cellsize=100e-6
thz=1e12
par='model=SphereArray resolution=2u simtime=100p resolution=2u radius=0u'

for wc in 2 4 8 16 32 64; do
	mpirun -np $NP  ../../scatter.py $par wirethick=1u wirecut=${wc}e-6
	../../effparam.py
done

sharedoptions='effparam/*.dat --paramname wirecut --figsizey 2 --xeval x/1e12 --ylim1 0 --parameval param*1e6'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
	--paramlabel 'wire cut $d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel 'wire cut $d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
	--paramlabel 'wire cut $d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag N' \
	--paramlabel 'wire cut $d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' \
	--paramlabel 'wire cut $d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_ni.pdf

