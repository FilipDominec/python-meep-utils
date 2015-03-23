#!/bin/bash

par='resolution=4u simtime=50p'

mpirun -np 1   python ../scatter.py $par wirethick=4u radius=0u  comment=OnlyWire  		             ;  ../effparam.py
mpirun -np 1   python ../scatter.py $par wirethick=0u radius=13u comment=OnlySphere  		         ;  ../effparam.py
mpirun -np 1   python ../scatter.py $par wirethick=4u radius=13u comment=SphereWireNIM  	         ;  ../effparam.py
mpirun -np 1   python ../scatter.py $par wirethick=4u radius=13u cellnumber=2 comment=SphereWireNIM  ;  ../effparam.py
mpirun -np 1   python ../scatter.py $par wirethick=4u radius=13u cellnumber=3 comment=SphereWireNIM  ;  ../effparam.py
