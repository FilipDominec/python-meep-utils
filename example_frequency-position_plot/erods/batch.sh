#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors

simtime=100e-12
resolution=2e-6
model=RodArray
cellnumber=3

for cellsize in 80e-6 90e-6 120e-6 
do
	## Preparation
	rm -r $model* effparam 2> /dev/null
	par="simtime=${simtime} cellsize=${cellsize} resolution=${resolution} model=$model orientation=E radius=10e-6 loss=.1 "

	## Generate the effective parameters
	mpirun -np $NP ../../scatter.py $par comment=nfile
	../../effparam.py
	 
	## Generate the field record over multiple cells
	mpirun -np $NP ../../scatter.py $par comment=fieldevolution cellnumber=${cellnumber} 
	../../effparam.py
	 
	## Plot the results
	../fxplot.py --nfile effparam/*nfile*dat  --fieldfile `cat last_simulation_name.dat`/FieldEvolution_at_*.h5 \
		--simtime ${simtime} --cellsize ${cellsize} --resolution ${resolution} --cellnumber ${cellnumber} \
		--freqlim1 0 --freqlim2 2e12 --numcontours 30 --fieldlim1 55% --fieldlim2 100% --Nconj yes --freqlim1 0.35e12 --freqlim2 1.25e12  \
		--decimatedata 3 --output erods_fxplot_cs${cellsize}.pdf 
done
