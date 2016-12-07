#!/bin/bash

## See the last line of ../metamaterial_models.py for the list of available models and their options
## XXX models = {'Slab':Slab, 'SphereArray':SphereArray, 'RodArray':RodArray}

## The computation procedure uses the 'cellsize' variable; make sure you set its correct physical value *before* 
## calling the compare_dispersion() procedure

if [ -z "$NP" ] ; then NP=2 ; fi             # number of processors
Ksampling=16
Kzones=2

compare_dispersion() {
	if [ -d cdh ]; then rm -r cdh/ ;  echo "Error: Another simulation has left a 'cdh/' directory; deleting it"; fi

    ## scan through the wave vector
    Kzs=`python -c "import numpy as np; print ' '.join(['%.4e' % x for x in np.linspace(0,$Kzones*np.pi/$cellsize,$Kzones*$Ksampling)])"`
	echo Will compute dispersion curves by finding all resonances at K in: $Kzs 
    for Kz in $Kzs
	do 
		echo "cs", $cellsize
        mpirun -np $NP ../cdh.py Kz=$Kz "$@" 
	done

    ../plot_cdh.py cdh/*dat --ylim2 $ylim2 ## (preview)

    ## compute the dispersion curves using the s-parameter method, (with the same parameters)
    mpirun -np $NP ../scatter.py  "$@" simtime=200p
    ../effparam.py

    ## repeat the plot, now comparing also to the curve retrieved above
    mv effparam/*dat NRef.dat
    ../plot_cdh.py cdh/*dat --ylim1 0.5 --ylim2 1.0 ## XXX

	## final cleanup
	lastname=`cat last_simulation_name.dat`
	mv cdh_ampli.png "CDH_${lastname}.png"
	mv cdh_ampli.pdf "CDH_${lastname}.pdf"
	mkdir "$lastname"
	mv cdh/ effparam/  NRef.dat "$lastname"
}


## === Rod array ===
ylim2=4 ## in THz
cellsize=100e-6
#for r in `seq 80 10 120`; do
#for r in 100; do
	#compare_dispersion ${par[@]} model=RodArray radius=${r}e-7
#done

## === Ordinary split-ring resonator === 
#cellsize=100e-6
#ylim2=4 ## in THz
#par=(resolution=4u simtime=30p cellsize=${cellsize})
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

## === Electro-magnetic symmetric split-ring resonator
#for icr in 6 8 10 18; do
for icr in 10 12 14 16 18; do
compare_dispersion ${par[@]} model=ESRRArray comment="emcSRR" cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
            insplitting=6e-6 incapacitorr=${icr}e-6 wirethick=0 radius=40e-6 srrthick=10e-6
done


#for icr in 12 14 16 18; do
#compare_dispersion ${par[@]} model=ESRRArray cbarthick=6e-6 splitting=6u  splitting2=6u capacitorr=5e-6 \
            #insplitting=6e-6 incapacitorr=${icr}e-6 wirethick=0 radius=30e-6 srrthick=10e-6
#done

#== Y-slits ==
#ylim2=4 ## in THz
#for cellsize in 50e-6 60e-6 75e-6; do
	#compare_dispersion ${par[@]} model=Fishnet slabthick=20u yholesize=inf xholesize=20u cellsizexy=100u cellsize=${cellsize}
#done
#mkdir YSlits; mv CDH* Fishnet* YSlits



## == Fishnets ==
## -- Testing our experimental fishnets made of steel -- 
#ylim2=1.5 ## in THz
#par=(resolution=4u simtime=100p cellsize=${cellsize}u )
#for cellsize in 100e-6; do # 50e-6 75e-6 100e-6 150e-6 200e-6 300e-6; do
	#compare_dispersion ${par[@]} model=Fishnet slabthick=20u xholesize=180u yholesize=200u cellsizexy=300u cellsize=${cellsize}
#done                                                                        
#mkdir RealisticFn_180x200_HR; mv CDH* Fishnet* RealisticFn_180x200_HR
#
#for cellsize in 100e-6; do                              # 50e-6 75e-6 100e-6 150e-6 200e-6 300e-6; do                             
	#compare_dispersion ${par[@]} model=Fishnet slabthick=20u xholesize=230u yholesize=255u cellsizexy=300u cellsize=${cellsize}
#done
#mkdir RealisticFn_230x255; mv CDH* Fishnet* RealisticFn_230x255

#for cellsize in 50e-6 100e-6 200e-6; do                             
	#compare_dispersion ${par[@]} model=Fishnet slabthick=20u xholesize=180u yholesize=200u cellsizexy=300u cellsize=${cellsize} slabcdist=30e-6
#done                                                                        
#mkdir DoubleRealisticFn_180x200; mv CDH* Fishnet* DoubleRealisticFn_180x200
#
#for cellsize in 50e-6 100e-6 200e-6; do                             
	#compare_dispersion ${par[@]} model=Fishnet slabthick=20u xholesize=230u yholesize=255u cellsizexy=300u cellsize=${cellsize} slabcdist=30e-6
#done
#mkdir DoubleRealisticFn_230x255; mv CDH* Fishnet* DoubleRealisticFn_230x255

#for cellsize in 50e-6 100e-6 200e-6; do                             
	#compare_dispersion ${par[@]} model=Fishnet slabthick=20u xholesize=200u yholesize=180u cellsizexy=300u cellsize=${cellsize}
#done
#mkdir RotatedRealisticFn_200x180; mv CDH* Fishnet* RotatedRealisticFn_200x180
#
#for cellsize in 50e-6 100e-6 200e-6; do                             
	#compare_dispersion ${par[@]} model=Fishnet slabthick=20u xholesize=255u yholesize=230u cellsizexy=300u cellsize=${cellsize}
#done
#mkdir RotatedRealisticFn_255x230; mv CDH* Fishnet* RotatedRealisticFn_255x230

## -- Navarro-Cia 2011 --
##experiment in microwave: 2.5x2.5 mm unit cell, stack periodicity 0.525 mm with 2.43 diel eps, hole radius 0.55 mm
#ylim2=3.5 ## in THz
#cellsize=52.5e-6
#par=(resolution=10u simtime=100p cellsize=${cellsize} )
#compare_dispersion ${par[@]} model=Fishnet slabthick=10u xholesize=110u yholesize=110u cornerradius=55u cellsizexy=250u resolution=5u 
#mkdir FnNavarroCia; mv CDH* Fishnet* FnNavarroCia

## -- Croenne and Wang 2010 --
#compare_dispersion ${par[@]} model=Fishnet slabthick=10u xholesize=110u yholesize=110u cornerradius=55u cellsizexy=250u resolution=5u  cellsize=52.5u
#mkdir CroenneWang; mv CDH* Fishnet* CroenneWang
