#!/bin/bash

## See the last line of ../metamaterial_models.py for the list of available models and their options
## XXX models = {'Slab':Slab, 'SphereArray':SphereArray, 'RodArray':RodArray}

if [ -z "$NP" ] ; then NP=2 ; fi             # number of processors

compare_dispersion() {
    ## scan through the wave vector
    for Kz in `seq 0  5000 62800` 
	do 
        mpirun -np $NP ../cdh.py Kz=$Kz simtime=30p "$@" 
	done

    ../plot_cdh.py cdh/*dat ## (preview)

    ## compute the dispersion curves using the s-parameter method, (with the same parameters)
    mpirun -np $NP ../scatter.py  "$@" simtime=200p padding=100u
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
par=(resolution=4u simtime=60p cellsize=$cellsize )
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

for icr in 12 14 16 18; do
compare_dispersion ${par[@]} model=ESRRArray cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
            insplitting=6e-6 incapacitorr=${icr}e-6 wirethick=0 radius=30e-6 srrthick=10e-6
done
