#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
cellsize=300e-6
thz=1e12
par="resolution=10u simtime=100p cellsize=$cellsize cellsizexy=$cellsize"

if [ -z "$skipsimulation" ]; then 
	echo
	mpirun -np $NP  ../../scatter.py $par comment="Simulation 180x200" model=Fishnet cornerradius=60u xholesize=180u yholesize=200u slabthick=10u
	../../effparam.py
	#mpirun -np $NP  ../../scatter.py $par comment="Simulation 230x255" model=Fishnet cornerradius=60u xholesize=230u yholesize=255u slabthick=10u
	#../../effparam.py
	#mpirun -np $NP  ../../scatter.py $par comment="Simulation 200x180" model=Fishnet cornerradius=60u xholesize=200u yholesize=180u slabthick=10u
	#../../effparam.py
	#mpirun -np $NP  ../../scatter.py $par comment="Simulation 255x230" model=Fishnet cornerradius=60u xholesize=255u yholesize=230u slabthick=10u
	#../../effparam.py
fi


sharedoptions=' --paramname comment --figsizey 2 --xeval x/1e12 --xlim2 1.5 --usetex yes'

../../plot_multiline.py $sharedoptions 28a_sec.dat 'effparam/Fishnet*180x200*_effparam.dat' --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '%s' --ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t_ap.pdf
../../plot_multiline.py $sharedoptions 28a_pri.dat 'effparam/Fishnet*200x180*_effparam.dat' --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '%s' --ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t_as.pdf
../../plot_multiline.py $sharedoptions 28c_sec.dat 'effparam/Fishnet*230x255*_effparam.dat' --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '%s' --ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t_cp.pdf
../../plot_multiline.py $sharedoptions 28c_pri.dat 'effparam/Fishnet*255x230*_effparam.dat' --xlabel "Frequency (THz)" --ycol '|t|' \
	--paramlabel '%s' --ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t_cs.pdf
