#!/bin/bash

## See the last line of ../metamaterial_models.py for the list of available models and their options
## XXX models = {'Slab':Slab, 'SphereArray':SphereArray, 'RodArray':RodArray}

if [ -z "$NP" ] ; then NP=2 ; fi             # number of processors


compare_dispersion() {
    ## scan through the wave vector
    for Kz in `seq 0  2000 60000`; do 
        mpirun -np $NP ../cdh.py Kz=$Kz $1; done
    ../plot_cdh.py cdh/*dat ## (preview)

    ## compute the dispersion curves using the s-parameter method, (with the same parameters)
    mpirun -np $NP ../scatter.py $1
    ../effparam.py

    ## repeat the plot, now comparing also to the curve retrieved above
    mv effparam/*dat NRef.dat
    ../plot_cdh.py cdh/*dat 
}


#compare_dispersion 'model=RodArray simtime=100p cellsize=100u'
#mv cdh_ampli.png RodArray_dispersion.png    # prevent the results from being overwritten
#mv cdh             cdh_RodArray            # archive the raw data

#compare_dispersion 'model=SphereArray comment=LoLoss wirethick=4u simtime=100p radius=30u cellsize=100u'
#mv cdh_ampli.png SphereArray_dispersion.png
#mv cdh             cdh_SphereArray

#compare_dispersion 'model=SphereArray resolution=6u comment=LoLoss wirethick=4u simtime=100p radius=40u cellsize=100u epsilon=12'
#mv cdh_ampli.png SphereArray_dispersion_r40_eps12.png
#mv cdh             cdh_SphereArray_r40_eps12
 
#compare_dispersion 'model=SRRArray simtime=100p'
#mv cdh_ampli.png    SRRArray_dispersion.png   
#mv cdh			    cdh_SRRArray             

compare_dispersion 'model=SRRArray simtime=100p splitting2=16u'
mpirun -np $NUMCPU ../scatter.py 'model=SRRArray simtime=100p'
mv cdh_ampli.png    SRRArrayDoubleSplit_dispersion.png   
mv cdh			    cdh_SRRArrayDoubleSplit             
