#!/usr/bin/env python
#-*- coding: utf-8 -*-

import numpy as np
## Harminv
def harminv(x, y, d=100, f=30, amplitude_prescaling=None):
    """
    suggested visualisation:
    hi = harminv(x,y)
    c = plt.scatter(hi['frequency'], hi['amplitude'], c=hi['phase'], s=np.abs(hi['quality'])/20 + 2, cmap=plt.cm.hsv, alpha=.3)
    """

    if amplitude_prescaling == None:
        amplitude_prescaling = 1000./np.max(np.abs(y))

    with open('/tmp/harminv_input.dat', 'w') as outfile: 
        outfile.write("#t[s]\t E(t)\n")
        np.savetxt(outfile, zip(x/2, np.real(y) * amplitude_prescaling), fmt="%.8e")
    import subprocess
    dt = x[1]-x[0]
    subprocess.Popen('harminv %f-%f -t %g   < /tmp/harminv_input.dat > /tmp/harminv_output.dat' % (0, 1./dt/10., dt*1.), shell=True, stdout=subprocess.PIPE).stdout.read().strip()
    try:
        (mf, md, mQ, mA, mp, merr) = np.loadtxt('/tmp/harminv_output.dat', usecols=list(range(6)), unpack=True, delimiter=', ', skiprows=1)
    except:
        print "\nWARNING: Harminv detected no resonances.\n\n"
        (mf, md, mQ, mA, mp, merr) = [np.array([]) for _ in range(6)]

    if type(mf) == np.float64:
        (mf, md, mQ, mA, mp, merr) = [np.array([val]) for val in (mf, md, mQ, mA, mp, merr)]

    return {'frequency':mf*2, 'decay':md * 2/np.pi, 'quality':mQ, 'amplitude':mA * 2 / amplitude_prescaling, 'phase':mp, 'error':merr}
