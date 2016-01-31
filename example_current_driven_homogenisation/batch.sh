#!/bin/bash

## See the last line of ../metamaterial_models.py for the list of available models and their options
## XXX models = {'Slab':Slab, 'SphereArray':SphereArray, 'RodArray':RodArray}

if [ -z "$NP" ] ; then NP=2 ; fi             # number of processors
Ksampling=10
Kzones=1

compare_dispersion() {
	if [ -d cdh ]; then echo "Error: Another simulation has left a 'cdh/' directory; delete it first"; exit; fi

    ## scan through the wave vector
    Kzs=`python -c "import numpy as np; print ' '.join(['%.4e' % x for x in np.linspace(0,$Kzones*np.pi/$cellsize,$Ksampling)])"`
	echo Will compute dispersion curves by finding all resonances at K in: $Kzs 
    for Kz in $Kzs
	do 
		echo "cs", $cellsize
        mpirun -np $NP ../cdh.py Kz=$Kz simtime=30p "$@" 
	done

    ../plot_cdh.py cdh/*dat ## (preview)

    ## compute the dispersion curves using the s-parameter method, (with the same parameters)
    mpirun -np $NP ../scatter.py  "$@" simtime=200p
    ../effparam.py

    ## repeat the plot, now comparing also to the curve retrieved above
    mv effparam/*dat NRef.dat
    ../plot_cdh.py cdh/*dat 

	## final cleanup
	lastname=`cat last_simulation_name.dat`
	mv cdh_ampli.png "CDH_${lastname}.png"
	mv cdh_ampli.pdf "CDH_${lastname}.pdf"
	mkdir "$lastname"
	mv cdh/ effparam/  NRef.dat "$lastname"
}

cellsize=100e-6
thz=1e12
par=(resolution=4u simtime=60p cellsize=${cellsize})
#compare_dispersion ${par[@]} model=SRRArray wirethick=0u radius=30u "comment=SRR only"
#compare_dispersion ${par[@]} model=SRRArray wirethick=4u radius=30u "comment=SRR with wire"
#compare_dispersion ${par[@]} model=SRRArray wirethick=0u radius=30u splitting2=16u "comment=symmetric SRR"
#compare_dispersion ${par[@]} model=SphereWire wirethick=0u radius=30u 'comment=TiO$_2$ sphere'
#compare_dispersion ${par[@]} model=SphereWire wirethick=4u radius=30u loss=.01 "comment=Lossless spheres with wires"
#compare_dispersion ${par[@]} model=SRRArray wirethick=4u radius=0u "comment=Wire only"
#compare_dispersion ${par[@]} model=SRRArray wirethick=4u radius=30u splitting2=16u "comment=symmetric SRR with wire"
#compare_dispersion ${par[@]} model=SRRArray wirethick=4u radius=30u splitting2=16u capacitorr=5e-6 "comment=symmetric SRR with pads"
#compare_dispersion ${par[@]} model=SRRArray resolution=4.000e-06 capacitorr=1.000e-05 splitting=4.000e-06 wirethick=4.000e-06 splitting2=4.000e-06 radius=3.000e-05 simtime=50p
#compare_dispersion ${par[@]} model=SRRArray resolution=4.000e-06 capacitorr=1.000e-05 splitting=4.000e-06 wirethick=0.000e-06 splitting2=4.000e-06 radius=3.000e-05 simtime=50p

#for icr in 6 8 10 18; do
#compare_dispersion ${par[@]} model=ESRRArray comment="emcSRR" cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
            #insplitting=6e-6 incapacitorr=${icr}e-6 wirethick=0 radius=40e-6 srrthick=10e-6
#done

#for r in `seq 80 10 120`; do
	#compare_dispersion ${par[@]} model=RodArray radius=${r}e-7
#done

#for icr in 12 14 16 18; do
#compare_dispersion ${par[@]} model=ESRRArray cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
            #insplitting=6e-6 incapacitorr=${icr}e-6 wirethick=0 radius=30e-6 srrthick=10e-6
#done


## == Fishnets ==
#thz=1e12
#par=(resolution=10u simtime=100p cellsize=${cellsize}u )
#for cellsize in 100e-6 150e-6 200e-6 300e-6; do
	#compare_dispersion ${par[@]} model=Fishnet slabcdist=0u slabthick=20u yholesize=200u xholesize=180u cellsizexy=300u resolution=10u cellsize=${cellsize}
#done
#for cellsize in 100e-6 150e-6 200e-6 300e-6; do
	#compare_dispersion ${par[@]} model=Fishnet slabcdist=0u slabthick=20u yholesize=255u xholesize=230u cellsizexy=300u resolution=10u cellsize=${cellsize}
#done

#for cellsize in 80e-6 70e-6  60e-6 50e-6; do
	#compare_dispersion ${par[@]} model=Fishnet slabcdist=30u slabthick=12u yholesize=50u xholesize=30u cellsizexy=100u cellsize=${cellsize}
#done

cellsize=80e-6
compare_dispersion ${par[@]} model=Fishnet slabcdist=0u slabthick=12u yholesize=30u xholesize=60u cellsizexy=100u cellsize=${cellsize}
compare_dispersion ${par[@]} model=Fishnet slabcdist=0u slabthick=12u yholesize=60u xholesize=60u cellsizexy=100u cellsize=${cellsize}

