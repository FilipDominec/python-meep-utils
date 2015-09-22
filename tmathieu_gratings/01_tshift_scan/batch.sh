#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
# par='model=TMathieu_Grating resolution=20n cellsize=100n'
par='model=TMathieu_Grating resolution=50n padding=10u ldist=15u tdist=10u rcore1=2u rcore2=2u'

for x in `seq 0 1 5`; do
	mpirun -np $NP   python ../../scatter.py $par simtime=500f tshift=${x}u
	../../effparam.py 
	echo 
done

../../plot_multiline.py TMathieu_Grating*.dat  \
		--xlabel 'Wavelength [$\upmu$m]' --xeval 1e6*c/x --xlim2=1.2   \
		--ycol '|t|' --yeval 'np.abs(y**2)'  --ylabel 'transmissivity $|s_{11}|^2$'   \
		--title 'Gold wire grids, $d_L$=15 $\mu$m, $d_T$=10 $\mu$m, $\rho_{1,2}$=2$\mu$m: a scan of the transverse shift $\Delta y$' \
        --paramname radius  --paramlabel '$\Delta y$ = %.0f $\mu$m'  --parameval 'param/1e-6'  \
		--usetex no  --output output_${PWD##*/}.png
