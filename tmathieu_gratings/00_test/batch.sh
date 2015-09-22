#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
# par='model=TMathieu_Grating resolution=20n cellsize=100n'
par='model=TMathieu_Grating resolution=50n padding=10u ldist=15u tdist=10u tshift=3u rcore1=2u rcore2=2u'

mpirun -np $NP   python ../../scatter.py $par simtime=300f
../../effparam.py 
