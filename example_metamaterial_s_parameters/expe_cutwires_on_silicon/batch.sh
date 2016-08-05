#!/bin/bash


for w in 30 40; do mpirun -np 2 ../../scatter.py model=WiresOnSi simtime=100p resolution=4u wirelength=${w}u; ../../effparam.py ; done

multiplot --ycol '|t|'  --ylabel 'Transmittance $|t|$'  --ylim2 1.0 --paramname comment --paramlabel %s --figsizey 3 \
	--xeval x/1e12 --xlim2 3.0  --xlabel "Frequency (THz)" \
	--output bousek.pdf     bousek_blue_1st_interface.dat WiresOnSi*dat 
