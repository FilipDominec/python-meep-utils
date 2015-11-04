#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
cellsize=100e-6
thz=1e12
par="model=SphereArray resolution=4u radius=30e-6 wirethick=0"

if [ -z "$skipsimulation" ]; then 
	for cellsize in `seq 60 10 200` `seq 220 20 400`; do
		mpirun -np $NP   python ../../scatter.py $par loss=1 simtime=150p  cellsize=${cellsize}u
		../../effparam.py
	done
fi

sharedoptions='effparam/*.dat --paramname cellsize --parameval param*1e6 --contours yes --numcontours 50 --colormap gist_earth --figsizex 4 --figsizey 3 --xeval x/1e12 --usetex yes'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
	--paramlabel 'a=%d $\upmu$m' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.png
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel 'a=%d $\upmu$m' \
   	--ylabel 'Transmittance $|t|$' --output ${PWD##*/}_t.png


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' --y2eval '0-y2' --ycol2 'imag N' \
	--paramlabel 'a=%d $\upmu$m' \
   	--ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.png 
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
	--paramlabel 'a=%d $\upmu$m' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.png  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag N' \
	--paramlabel 'a=%d $\upmu$m' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.png


#../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' --y2eval '0-y2' --ycol2 'imag eps' \
	#--paramlabel 'a=%d $\upmu$m' \
   	#--ylabel 'Permittivity $\varepsilon_{\text{eff}}$' --output ${PWD##*/}_eps.png
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' \
	--paramlabel 'a=%d $\upmu$m' --ylim1 -5 --ylim2 10 \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_epsr.png
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag eps' \
	--paramlabel 'a=%d $\upmu$m' \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_epsi.png
#
#
#../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' --y2eval '0-y2' --ycol2 'imag mu' \
	#--paramlabel 'a=%d $\upmu$m' \
   	#--ylabel 'Permeability $\mu_{\text{eff}}$' --output ${PWD##*/}_mu.png
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' \
	--paramlabel 'a=%d $\upmu$m' \
   	--ylabel 'Permeability $\mu_{\text{eff}}^{\prime}$' --output ${PWD##*/}_mur.png
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag mu' \
	--paramlabel 'a=%d $\upmu$m' \
   	--ylabel 'Permeability $\mu_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_mui.png
