#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi             # number of processors
if [ -z $ext ] ; then ext=pdf ; fi             # number of processors
thz=1e12
cellsize=100e-6
staticpar=(model=Fishnet resolution=2u simtime=150p cellsize=100u cellsizexy=100u slabcdist=0u slabthick=10e-6 ) 

if [ -z "$skipsimulation" ]; then 
    for    xhs in `seq 2 8 45` `seq 6 8 100`; do
        np ../../scatter.py "${staticpar[@]}" yholesize=inf xholesize=${xhs}u
        ../../effparam.py --numstabcoef .97
    done
fi

plotoptions=(effparam/*.dat  --xlabel     'Frequency (THz)'       --xeval x/1e12 \
        --paramname xholesize  --paramlabel 'Slit width $d_x$' --parameval param*1e6 \
        --contours yes  --numcontours 25  --colormap gist_earth  --figsizex 4  --figsizey 3)

../../plot_multiline.py "${plotoptions[@]}" --ycol '|r|'    --ylabel 'Reflectance   $|r|$'  --ylim1 0 --ylim2 1 \
        --output ${PWD##*/}_r.$ext 
../../plot_multiline.py "${plotoptions[@]}" --ycol '|t|'    --ylabel 'Transmittance $|t|$'  --ylim1 0 --ylim2 1 \
        --output ${PWD##*/}_t.$ext 
../../plot_multiline.py "${plotoptions[@]}" --ycol 'real N' --ylabel 'Refractive index $N_{\text{eff}}^{\prime}$' \
        --output ${PWD##*/}_nr.$ext
../../plot_multiline.py "${plotoptions[@]}" --ycol 'imag N' --ylabel 'Refractive index $N_{\text{eff}}^{\prime\prime}$' --yeval '0-y' \
        --output ${PWD##*/}_ni.$ext  --ylim1 -1 --ylim2 0

#../../plot_multiline.py "${sharedoptions[@]}" --ycol 'real eps' --ycol2 'imag eps' \
       #--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime}$' --output ${PWD##*/}_eps.pdf  --y2eval '0-y2' --ylim1 -12 --ylim2 3 
#../../plot_multiline.py "${sharedoptions[@]}" --yeval '0-y' --ycol 'imag eps' --ylim1 -12 --ylim2 3 \
       #--ylabel 'Effective permittivity $\varepsilon_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_epsi.pdf
#
#../../plot_multiline.py "${sharedoptions[@]}" --ycol 'real mu' --y2eval '0-y2' --ycol2 'imag mu' --ylim1 -5 --ylim2 5 \
    #-ylabel 'Permeability $\mu_{\text{eff}}$' --output ${PWD##*/}_mu.pdf
#../../plot_multiline.py "${sharedoptions[@]}" --ycol 'real mu' --ylim1 -5 --ylim2 5 \
    #--ylabel 'Permeability $\mu_{\text{eff}}^{\prime}$' --output ${PWD##*/}_mur.pdf
#../../plot_multiline.py "${sharedoptions[@]}" --yeval '0-y' --ycol 'imag mu' --ylim1 -5 --ylim2 5 \
    #--ylabel 'Permeability $\mu_{\text{eff}}^{\prime\prime}$' --output ${PWD##*/}_mui.pdf
