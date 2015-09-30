#!/bin/bash
if [ -z $NP ] ; then NP=2 ; fi			 # number of processors
mpirun -np $np  ../scatter.py model=HalfSpace blend=10u comment=metal ; ../effparam.py
