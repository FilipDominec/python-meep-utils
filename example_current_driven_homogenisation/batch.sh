#!/bin/bash

for Kz in `seq 0    5000 150000`; do ../cdh.py Kz=$Kz simtime=30p; done
../plot_cdh.py S*dat ## preview
for Kz in `seq 2500 5000 150000`; do ../cdh.py Kz=$Kz simtime=30p; done
../plot_cdh.py S*dat

