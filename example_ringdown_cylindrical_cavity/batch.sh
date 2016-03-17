#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors

## Settings
radius=33.774e-3
height=122.36e-3


## perform the actual electrodynamic computation
mpirun -np $NP ../cylindrical_cavity.py	radius=$radius height=$height	

## save `annotate.txt` containing the resonance mode frequencies
./cylindrical_cavity_modes_annotation.py --radius $radius --height $height					

## analyze the time-domain waveform
../ringdown_analysis.py HollowCyl_timedomain.dat			
mv annotate.txt annotate_cylinder.txt
mv oscillator_spectrum.png Cavity_oscillator_spectrum_height_${height}_radius_${radius}.png


## Optionally look how Harminv can resolve experimental terahertz resonance peaks of water vapour
#cp annotate_water.txt annotate.txt
#../ringdown_analysis.py input_tdts.dat
#mv oscillator_spectrum.png terahertz_H2Ovapour_spectrum.png
