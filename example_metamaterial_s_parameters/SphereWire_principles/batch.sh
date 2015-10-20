#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
cellsize=100e-6
thz=1e12
par='model=SphereArray resolution=4u simtime=100p'

#mpirun -np $NP  ../../scatter.py $par wirethick=4u radius=13u comment='TiO$_2$ spheres with wires'
#../../effparam.py
#mpirun -np $NP  ../../scatter.py $par wirethick=4u radius=13u comment='Lossless spheres with wires' loss=.01
#../../effparam.py
#mpirun -np $NP  ../../scatter.py $par wirethick=4u radius=0u  comment='Wires only'
#../../effparam.py
#mpirun -np $NP  ../../scatter.py $par wirethick=0u radius=13u comment='TiO$_2$ Spheres only'
#../../effparam.py

sharedoptions='effparam/*.dat --paramname comment --figsizey 2 --xeval x/1e12 --ylim1 0 --xlim2 1.5'

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf # --color RdYlBu

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf  #--color RdYlBu_r

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^\prime$' --output ${PWD##*/}_nr.pdf  \
    --overlayplot "-c/2/$cellsize/x/$thz, c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz"  

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_ni.pdf #--color PiYG_r

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real eps' \
	--paramlabel 'wire cut $d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_epsr.pdf

../../plot_multiline.py $sharedoptions --xlabel "Frequency (THz)" --ycol 'real mu' \
	--paramlabel 'wire cut $d_c = %.0f$ $\upmu$m' \
   	--ylabel 'Permeability $\mu_{\text{eff}}^{\prime}$' --output ${PWD##*/}_mur.pdf


