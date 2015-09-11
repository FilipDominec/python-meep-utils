#!/bin/bash
if [ -z $NP ] ; then NP=1 ; fi			 # number of processors
par='model=SphereArray resolution=4u simtime=100p'

for loss in 0.001 0.01 0.03 0.1 0.3 1; do 
    mpirun -np $np   python ../../scatter.py $par loss=$loss; ../../effparam.py
done
../../plot_multiline.py effparam/SphereArray_simtime\=1.000e-10_*.dat --paramname loss --paramlabel 'Losses percentage %d' --paramunit .01  --ycol 'real N' --xunit 1e12
