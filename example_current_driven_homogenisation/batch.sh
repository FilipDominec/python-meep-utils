#!/bin/bash

PARAM='comment=LoLoss wirethick=4u simtime=100p'

for Kz in `seq 0  5000 60000`; do ../cdh.py Kz=$Kz $PARAM; done
../plot_cdh.py cdh/*dat ## preview


../scatter.py $PARAM
../effparam.py
mv effparam/*dat NRef.dat

../plot_cdh.py cdh/*dat ## now plot with the curve retrieved with s-parameters
