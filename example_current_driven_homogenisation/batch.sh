#!/bin/bash

## See the last line of ../metamaterial_models.py for the list of available models and their options
## XXX models = {'Slab':Slab, 'SphereArray':SphereArray, 'RodArray':RodArray}

if [ -z "$NUMCPU" ]; then 
	NUMCPU=1; 
	echo "Note:Defaulting to one processor, use e.g. 'export NUMCPU=4' to use multiprocessing"; fi

compare_dispersion() {
	## scan through the wave vector
	for Kz in `seq 0  5000 60000`; do 
		mpirun -np $NUMCPU ../cdh.py Kz=$Kz $1; done
	../plot_cdh.py cdh/*dat ## (preview)

	## compute the dispersion curves using the s-parameter method, (with the same parameters)
	mpirun -np $NUMCPU ../scatter.py $1
	../effparam.py

	## repeat the plot, now comparing also to the curve retrieved above
	mv effparam/*dat NRef.dat
	../plot_cdh.py cdh/*dat 
}

compare_dispersion 'model=SphereArray comment=LoLoss wirethick=4u simtime=100p'
mv cdh SphereArray_results		# prevent the results from being overwritten

compare_dispersion 'model=RodArray comment=LoLoss wirethick=4u simtime=100p'
mv cdh RodArray_results			# prevent the results from being overwritten
