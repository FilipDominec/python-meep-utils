#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
par='model=SphereArray resolution=4u simtime=100p'

for loss in 0.003 0.01 0.03 0.1 0.3 1; do 
    mpirun -np $np   python ../../scatter.py $par loss=$loss; ../../effparam.py
done

../effparam.py

../../plot_multiline.py effparam/SphereArray_simtime\=1.000e-10_*.dat --paramname loss --paramlabel 'Losses percentage %d' --paramunit .01  --ycol 'real N' --xlabel 'Frequency [THz]' --xeval x/1e12

../../plot_multiline.py effparam/RodArray_effparam.dat    --paramname radius --paramlabel none \
        --xlabel "Frequency (THz)" --xeval x/1e12  \
        --ycol '|r|' --ylabel 'Reflectance   $|r|$' --figsizey 2 --output ${PWD##*/}_r.pdf --color RdYlBu

../../plot_multiline.py effparam/RodArray_effparam.dat    --paramname radius --paramlabel none  \
        --xlabel "Frequency (THz)" --xeval x/1e12  \
        --ycol '|t|' --ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf --color RdYlBu_r

../../plot_multiline.py effparam/RodArray_effparam.dat    --paramname radius --paramlabel none  \
        --xlabel "Frequency (THz)" --xeval x/1e12  \
        --ycol 'real N' --ylim2 5.  --ylabel 'Effective index $N^\prime$' --figsizey 2 --output ${PWD##*/}_n.pdf --color PiYG_r \
		--overlayplot "2.998e8/4./50e-6/(x*1e12)"

