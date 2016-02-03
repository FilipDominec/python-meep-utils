#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
cellsize=100e-6
simtime=200e-12
resolution=4e-6
par="simtime=${simtime} cellsize=${cellsize} resolution=${resolution} model=RodArray orientation=E radius=10e-6 loss=.1 "

## Generate the effective parameters
mpirun -np $NP ../../scatter.py $par comment=nfile
../../effparam.py
 
## Generate the field record over multiple cells
cellnumber=3
mpirun -np $NP ../../scatter.py $par comment=fieldfile cellnumber=${cellnumber} 
../../effparam.py
 
## Plot the results
../fxplot.py --nfile effparam/*nfile*dat  --fieldfile `cat last_simulation_name.dat`/FieldEvolution_at_*.h5 \
    --simtime ${simtime} --cellsize ${cellsize} --resolution ${resolution} --cellnumber ${cellnumber} \
   	--freqlim1 0 --freqlim2 2e12 --fieldlim1 20% --fieldlim2 95% --output erods_f-x_plot.pdf
