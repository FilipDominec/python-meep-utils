#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
if [ -z $ext ] ; then ext=png ; fi             # number of processors
thz=1e12
cellsize=100e-6
staticpar=(model=Fishnet resolution=1u simtime=50p cellsize=100u cellsizexy=100u slabcdist=0u xholesize=20u)

if [ -z "$skipsimulation" ]; then 
    for K in `seq 0 500000 5000000` 50000 200000 `seq 6000000 1000000 12000000`; do
        #np ../../scatter.py "${staticpar[@]}" 
        mpirun -np $NP  ../scatter.py model=HalfSpace simtime=500f resolution=15n blend=0u padding=1u comment=metal-fieldevolution Kx=${K} 
        ../effparam.py
    done
fi

plotoptions=(Half*.dat  --xlabel     'Frequency'   \
        --paramname Kx  --paramlabel 'Angle $a$' --parameval "np.arcsin(param*c/2/np.pi/x)/np.pi*180" \
        --contours yes  --numcontours 25  --colormap gist_earth  --figsizex 4  --figsizey 3)
  


../plot_multiline.py "${plotoptions[@]}" --ycol '|r|'    --ylabel 'Reflectance   $|r|$'  --ylim1 0 --ylim2 1 \
        --output ${PWD##*/}_r.$ext 
#../../plot_multiline.py "${plotoptions[@]}" --ycol '|t|'    --ylabel 'Transmittance $|t|$'  --ylim1 0 --ylim2 1 \
        #--output ${PWD##*/}_t.$ext 
