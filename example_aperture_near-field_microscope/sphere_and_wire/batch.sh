#!/bin/bash

np=2
otherparams='simtime=30p'

mpirun -np $np ../../nearfield.py radius=-1 wireth=-1 comment=Ref  $otherparams 
mv Apert*comment=Ref*.dat ref.dat

mpirun -np $np ../../nearfield.py radius=10u wireth=-1 $otherparams 
