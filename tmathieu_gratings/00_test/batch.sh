#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
# par='model=TMathieu_Grating resolution=20n cellsize=100n'
par='model=TMathieu_Grating resolution=50n padding=10u ldist=15u tdist=10u tshift=3u rcore1=2u rcore2=2u'

#for radius in `seq 0 5 50`; do
	mpirun -np $NP   python ../../scatter.py $par simtime=300f
	../../effparam.py 
	#echo 
#done

#../plot_multiline.py PlasmonicSpheres_*.dat  \
		#--xlabel 'Wavelength [$\upmu$m]' --xeval 1e6*c/x --xlim2=1.2   \
		#--ycol '|r|' --yeval 'np.abs(y**2)'  --ylabel 'reflectivity $|s_{11}|^2$'   \
		#--title 'Gold spheres in 40nm teflon, 100$\times$100nm cell: a scan of the radius $\rho$' \
        #--paramname radius  --paramlabel '$\rho$ = %.0f nm'  --parameval 'param/1e-9'  \
		#--usetex yes  --output output_${PWD##*/}.png
