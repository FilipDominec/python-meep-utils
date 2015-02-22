#!/bin/bash
#mpirun -np 1   python ../scatter.py resolution=4u simtime=50p padding=00u


for ff in `seq 90 100`; do  
    mpirun -np 1   python ../scatter.py resolution=4u simtime=50p padding=00u frequency=${ff}0e9
done


