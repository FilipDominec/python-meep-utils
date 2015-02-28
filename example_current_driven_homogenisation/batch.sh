#!/bin/bash

for Kz in `seq 0    5000 60000`; do ../cdh.py Kz=$Kz simtime=30p; done
../plot_cdh.py cdh/*dat ## preview


../scatter.py comment=LossLess
../effparam.py
mv effparam/*dat NRef.dat

../plot_cdh.py cdh/*dat ## now plot with the curve retrieved with s-parameters
