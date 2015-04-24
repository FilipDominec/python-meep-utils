#!/bin/bash

np=2
otherparams='simtime=50p'

rm -f ref.dat
mpirun -np $np ../../nearfield.py radius=-1 wireth=-1 comment=Ref  $otherparams 
mv Apert*comment=Ref*.dat ref.dat

mpirun -np $np ../../nearfield.py radius=10u wireth=-1 $otherparams  comment=Sphere
mpirun -np $np ../../nearfield.py radius=-1  wireth=4u $otherparams  comment=Wire 
mpirun -np $np ../../nearfield.py radius=10u wireth=4u $otherparams  comment='Sphere with a wire'

../../plot_multiline.py *NORMALIZED* --paramname comment
