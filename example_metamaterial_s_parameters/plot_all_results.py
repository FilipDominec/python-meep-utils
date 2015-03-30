#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, hbar, pi

## Use LaTeX
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', size=12)
matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{yfonts}, \usepackage{txfonts}') #, \usepackage{palatino} \usepackage{lmodern}, 
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})  ## select fonts

## Start figure + subplot 
fig = plt.figure(figsize=(10,10))
#fig.subplots_adjust(left=.05, bottom=.05, right=.99, top=.99, wspace=.05, hspace=.05) ## (for interactive mode)


## Sort arguments by a _numerical_ value in their parameter, keep the color order
filenames = sys.argv[1:]
paramname = 'wirethick'
def get_param(filename):             ## Load header to the 'parameters' dictionary
    parameters = {}
    with open(filename) as datafile:
        for line in datafile:
            if (line[0:1] in '0123456789') or ('column' in line.lower()): break    # end of parameter list
            key, value = line.replace(',', ' ').split()[-2:]
            try: value = float(value) ## Try to convert to float, if possible
            except: pass                ## otherwise keep as string
            parameters[key] = value
    return parameters
params  = [get_param(n)[paramname] for n in filenames]
datasets = zip(params, filenames)                               ## sort the files by the parameter
datasets.sort()
colors = matplotlib.cm.gist_earth(np.linspace(0.1,0.9,len(filenames)+1)[:-1]) ## add the colors to sorted files
datasets = zip(colors, *zip(*datasets))


freq_unit = 1e12

ax = plt.subplot(311, axisbg='w')
for color, param, filename in datasets:
    (x, r) = np.loadtxt(filename, usecols=[0,1], unpack=True)
    plt.plot(x/freq_unit, r,color=color, label='$p = %.3g$' % (param))
plt.ylabel(u"y"); 
plt.grid()
plt.legend(prop={'size':12}, loc='upper left') #.draw_frame(False)


ax = plt.subplot(312, axisbg='w')
plt.axvspan(1e12/freq_unit, 1.5e12/freq_unit, color='b', alpha=.1); 
for color, param, filename in datasets:
    (x, t) = np.loadtxt(filename, usecols=[0,3], unpack=True)
    plt.plot(x/freq_unit, t,color=color, label='$p = %.3g$' % (param))
plt.ylim((1e-3, 1.1e0))
plt.yscale('log')
plt.ylabel(u"y"); 
plt.grid()
plt.legend(prop={'size':12}, loc='upper left') #.draw_frame(False)

ax = plt.subplot(313, axisbg='w')
for color, param, filename in datasets:
    (x, r, t) = np.loadtxt(filename, usecols=[0,1,3], unpack=True)
    plt.plot(x/freq_unit, 1-r**2-t**2,color=color, label='ko\\v ci\\v cka $\int\lambda_{\Gamma} x = %.3g$' % (param))
plt.xlabel(u"frequency [THz]"); 
plt.ylabel(u"y"); 
plt.grid()
plt.legend(prop={'size':12}, loc='upper left') #.draw_frame(False)


## ==== Outputting ====
plt.savefig("output.png", bbox_inches='tight')
plt.savefig("output.pdf", bbox_inches='tight')
