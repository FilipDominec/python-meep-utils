for x in `seq 5 5 50` ; do np ../../scatter.py model=PlasmonicDimers radius=${x}e-10 simtime=200f resolution=.5n comment=noosc-singlesphere ; done
