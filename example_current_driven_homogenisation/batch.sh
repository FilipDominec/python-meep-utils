#!/bin/bash

## See the last line of ../metamaterial_models.py for the list of available models and their options
## XXX models = {'Slab':Slab, 'SphereArray':SphereArray, 'RodArray':RodArray}

if [ -z "$NP" ] ; then NP=2 ; fi             # number of processors

compare_dispersion() {
    ## scan through the wave vector
    for Kz in `seq 0  2000 62800` 
	do 
        mpirun -np $NP ../cdh.py Kz=$Kz "$@" 
	done

    ../plot_cdh.py cdh/*dat ## (preview)

    ## compute the dispersion curves using the s-parameter method, (with the same parameters)
    mpirun -np $NP ../scatter.py "$@"
    ../effparam.py

    ## repeat the plot, now comparing also to the curve retrieved above
    mv effparam/*dat NRef.dat
    ../plot_cdh.py cdh/*dat 

	## final cleanup
	lastname=`cat last_simulation_name.dat`
	mv cdh_ampli.png    "CDH${lastname}.png"
	rm -r cdh/ effparam/  NRef.dat
}

cellsize=100e-6
thz=1e12
par=(resolution=4u simtime=30p cellsize=$cellsize )
compare_dispersion ${par[@]} model=ESRRArray comment="emcSRR" cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
            insplitting=6e-6 incapacitorr=10e-6 wirethick=0 radius=40e-6 srrthick=10e-6
compare_dispersion ${par[@]} model=ESRRArray comment="emcSRR" cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
            insplitting=6e-6 incapacitorr=12e-6 wirethick=0 radius=40e-6 srrthick=10e-6
compare_dispersion ${par[@]} model=ESRRArray comment="emcSRR" cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
            insplitting=6e-6 incapacitorr=8e-6 wirethick=0 radius=40e-6 srrthick=10e-6
compare_dispersion ${par[@]} model=SRRArray wirethick=0u radius=30u "comment=SRR only"
compare_dispersion ${par[@]} model=SRRArray wirethick=4u radius=30u "comment=SRR with wire"
compare_dispersion ${par[@]} model=SRRArray wirethick=0u radius=30u splitting2=16u "comment=symmetric SRR"
compare_dispersion ${par[@]} model=SphereWire wirethick=0u radius=30u 'comment=TiO$_2$ sphere'
compare_dispersion ${par[@]} model=SphereWire wirethick=4u radius=30u loss=.01 "comment=Lossless spheres with wires"
#compare_dispersion ${par[@]} model=SRRArray wirethick=4u radius=0u "comment=Wire only"

