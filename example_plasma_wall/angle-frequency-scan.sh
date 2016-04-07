#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
if [ -z $ext ] ; then ext=png ; fi             # number of processors

staticpar=(model=HalfSpace simtime=250f resolution=15n padding=.5u epsilon=2)

if [ -z "$skipsimulation" ]; then 
    #for K in 0.1 0.3 `seq 0 1 9` `seq 10 2 20`; do    ## transverse wavenumber in 1/um
    for K in 0.1 0.3 `seq 0 1 9` `seq 10 2 20` `seq .5 1 9 | tr , .` `seq 11 2 20`; do    ## transverse wavenumber in 1/um
        mpirun -np $NP  ../../scatter.py "${staticpar[@]}" blend=${blend} comment=${comment} epsilon=2 $Kcomponent=${K}e6
    done
fi

## Resolve the angle of incidence from the magnitude of the wavevector and its transverse component and plot
plotoptions=(Half*.dat  --xlabel     'Frequency'   \
        --paramname Kx  --paramlabel 'Angle $a$' --parameval "np.arcsin(param*c/2/np.pi/x)/np.pi*180" \
        --contours yes  --numcontours 25  --colormap gist_earth  --figsizex 4  --figsizey 3 --interp_aspect .5)

../../plot_multiline.py "${plotoptions[@]}" --ycol '|r|'    --ylabel 'Reflectance   $|r|$'  --ylim1 0 --ylim2 1 \
        --output ../${PWD##*/}_r.$ext 

## Clean up
if [ -z "$skipsimulation" ]; then 
    rm Half*/
fi
