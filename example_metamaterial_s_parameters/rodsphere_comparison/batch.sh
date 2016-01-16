#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
thz=1e12
cellsize=100e-6
par="cellsize=$cellsize simtime=150p"
if [ -z "$skipsimulation" ]; then 
	#mpirun -np $NP ../../scatter.py $par model=RodArray resolution=1.5u orientation=H radius=19e-6 comment='Rods $||\;\mathbf{H}$'
	#../../effparam.py
	mpirun -np $NP ../../scatter.py $par model=SphereWire resolution=3u radius=30e-6 comment='Spheres'
	../../effparam.py
fi

sharedoptions="effparam/*.dat --paramname comment  --parameval param*1e6 --figsizey 2 --xeval x/1e12"
../../plot_multiline.py $sharedoptions --paramlabel '%s' --xlabel "Frequency (THz)" --ycol '|r|' \
   	--ylabel 'Reflectance   $|r|$' --output ${PWD##*/}_r.pdf
../../plot_multiline.py $sharedoptions --paramlabel '%s' --xlabel "Frequency (THz)" --ycol '|t|' \
   	--ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf 


../../plot_multiline.py $sharedoptions --paramlabel '%s' --xlabel "Frequency (THz)" --ycol 'real N' --y2eval '0-y2' --ycol2 'imag N' \
   	--ylabel 'Refractive index $N_{\text{eff}}$' --output ${PWD##*/}_n.pdf  \
    --overlayplot "c/2/$cellsize/x/$thz,2*c/2/$cellsize/x/$thz,3*c/2/$cellsize/x/$thz,4*c/2/$cellsize/x/$thz" --ylim1 -2 --ylim2 2.5  

#../../plot_multiline.py $sharedoptions --paramlabel '%s' --xlabel "Frequency (THz)" --ycol 'real eps' --y2eval '0-y2' --ycol2 'imag eps' --ylim1 -12 --ylim2 3  \
   	#--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_eps.pdf \
    #--overlayplot "1-.68**2/x**2"  
#../../plot_multiline.py $sharedoptions --paramlabel '%s' --xlabel "Frequency (THz)" --ycol 'real eps' --ylim1 -12 --ylim2 3  \
   	#--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_epsr.pdf \
    #--overlayplot "1-.68**2/x**2"  
#../../plot_multiline.py $sharedoptions --paramlabel '%s' --xlabel "Frequency (THz)" --yeval '0-y' --ycol 'imag eps' --ylim1 -12 --ylim2 3  \
   	#--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_epsi.pdf \
#
