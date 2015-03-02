#!/bin/bash
#mpirun -np 1 ../plasmons.py comment=Asym metalthick=100n resolution=50n size_x=8e-6 apdistx=6e-6 simtime=70f
mpirun -np 1 ../plasmons.py comment=Asym metalthick=100n resolution=50n size_x=8e-6 apdistx=6e-6 simtime=20f
