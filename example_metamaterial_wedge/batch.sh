#!/bin/bash

## example:
## ./batch.sh radius=10u resolution=4u simtime=30p

../scatter.py model=RodArray $@
../effparam.py
last_simulation_name=`cat last_simulation_name.dat`
rm -f effparam.dat; ln -s effparam/${last_simulation_name}_effparam.dat effparam.dat

../wedge_refraction.py $@
last_simulation_name=`cat last_simulation_name.dat`
../scripts_postpro/plot_FA.py $last_simulation_name/wedge_at_x0.000e+00_z*.h5
