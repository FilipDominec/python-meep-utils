#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
## Harminv
def harminv(x,y, d=100, f=30):
    """
    suggested visualisation:
    hi = harminv(x,y)
    c = plt.scatter(hi['frequency'], hi['amplitude'], c=hi['phase'], s=np.abs(hi['quality'])/20 + 2, cmap=plt.cm.hsv, alpha=.3)
    """


    with open('hitest.dat', 'w') as outfile: 
        outfile.write("#t[s]\t E(t)\n")
        np.savetxt(outfile, zip(x/2, y), fmt="%.8e")
    import subprocess
    dt = x[1]-x[0]
    subprocess.Popen('harminv 0-1000 -t %g -d %g -f %g < hitest.dat > hiout.dat' % (dt, d, f), shell=True, stdout=subprocess.PIPE).stdout.read().strip()
    (mf, md, mQ, mA, mp, merr) = np.loadtxt('hiout.dat', usecols=list(range(6)), unpack=True, delimiter=', ', skiprows=1)
    return {'frequency':mf*2, 'decay':md, 'quality':mQ, 'amplitude':mA, 'phase':mp, 'error':merr}
