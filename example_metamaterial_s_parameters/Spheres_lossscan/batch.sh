#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
cellsize=100e-6
thz=1e12
par="model=SphereArray resolution=4u radius=30e-6 wirethick=0 cellsize=$cellsize"

mpirun -np $NP   python ../../scatter.py $par loss=1 comment='100' simtime=150p
../../effparam.py
mpirun -np $NP   python ../../scatter.py $par loss=0.3	comment='10' simtime=150p
../../effparam.py
mpirun -np $NP   python ../../scatter.py $par loss=0.1 comment='1' simtime=300p
../../effparam.py

sharedoptions='effparam/*.dat --paramname comment --figsizey 2 --xeval x/1e12 --ylim1 0 --xlim2 1.5'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' --ycol2 'imag N' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag N' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' --ycol2 'imag eps' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}$' --output ${PWD##*/}_eps.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_epsr.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag eps' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_epsi.pdf


../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' --ycol2 'imag mu' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Permeability $\mu_{\text{eff}}$' --output ${PWD##*/}_mu.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Permeability $\mu_{\text{eff}}^{\prime}$' --output ${PWD##*/}_mur.pdf
../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag mu' \
	--paramlabel '%.1f\%% losses' \
   	--ylabel 'Permeability $\mu_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_mui.pdf
