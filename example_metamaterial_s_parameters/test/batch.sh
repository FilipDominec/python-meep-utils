#!/bin/batch
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
model=SphereArray
cellsize=50e-6
mpirun -np $NP ../scatter.py model=SphereArray resolution=4u simtime=100p wirethick=10u cellsize=$cellsize padding=20e-6 radius=13e-6
../../effparam.py 
../../plot_multiline.py RodArray.dat    --paramname radius --paramlabel '' --ycol '|r|' --contours no  \
                    --yeval 'np.abs(y)' --ylim1 0 --ylim2 1 \
                    --xlabel 'Normalized grating period $\Lambda/\lambda$' --xeval x*1e-9/0.8e-6 --xlim1=0    \
                    --paramlabel 'Reflectance $\muup$m)' --parameval param*1e-9/1e-6  \
                    --title "Reflectance $|r|^2$ of ${comment%%-narrowfreq}$\\perp$ at $\\lambda$=800 nm)" --usetex yes --colormap 'Paired/2' \
                    --output ${PWD##*/}_r.png

../effparam.py

../../plot_multiline.py effparam/RodArray_effparam.dat    --paramname radius --paramlabel none \
        --xlabel "Frequency (THz)" --xeval x/1e12  \
        --ycol '|r|' --ylabel 'Reflectance   $|r|$' --figsizey 2 --output ${PWD##*/}_r.pdf --color RdYlBu

../../plot_multiline.py effparam/RodArray_effparam.dat    --paramname radius --paramlabel none  \
        --xlabel "Frequency (THz)" --xeval x/1e12  \
        --ycol '|t|' --ylabel 'Transmittance $|t|$' --figsizey 2 --output ${PWD##*/}_t.pdf --color RdYlBu_r

../../plot_multiline.py effparam/RodArray_effparam.dat    --paramname radius --paramlabel none  \
        --xlabel "Frequency (THz)" --xeval x/1e12  \
        --ycol 'real N' --ylim2 5.  --ylabel 'Effective index $N^\prime$' --figsizey 2 --output ${PWD##*/}_n.pdf --color PiYG_r \
		--overlayplot "2.998e8/4./50e-6/(x*1e12)"


