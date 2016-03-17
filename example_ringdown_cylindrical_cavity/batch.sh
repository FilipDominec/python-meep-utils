#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors

echo == Resonances in cylindrical cavity ==
mpirun -np $NP ../cylindrical_cavity.py									## performs the actual electrodynamic computation
./cylindrical_cavity_modes_annotation.py					## saves `annotate.txt` containing the resonance mode frequencies
../ringdown_analysis.py HollowCyl_timedomain.dat			## analyzes the time-domain waveform
mv annotate.txt annotate_cylinder.txt
mv oscillator_spectrum.png Cavity_oscillator_spectrum.png

echo 
echo == Experimental resonance peaks of water vapour ==
cp annotate_water.txt annotate.txt
../ringdown_analysis.py input_tdts.dat
mv oscillator_spectrum.png terahertz_H2Ovapour_spectrum.png
