#!/usr/bin/env python
#-*- coding: utf-8 -*-

dataset = 'all'
#dataset = 'ex.r'
#dataset = 'ex.i'

## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, hbar, pi

## Use LaTeX
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', size=8)
matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{yfonts}, \usepackage{txfonts}, \usepackage{lmodern},')
## LaTeX options
#matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})  ## select fonts
#matplotlib.rc('text.latex',unicode=True)   ## use accented letters


## === 2. Load data ===

## Load field arrays from a HDF5 file
import h5py
h5file = h5py.File('snapshot_Ex_t0.000000e+00.h5', "r")
print "Found datasets:", h5file.keys()
if dataset != 'all':
    data = np.array(h5file[dataset])[0,:,:,0]
else:
    data = np.array(h5file['ex.r'])[0,:,:,0] + 1j*np.array(h5file['ex.i'])[0,:,:,0]
print "Loaded dataset with shape:", data.shape
z = data
x = np.linspace(0, 1, z.shape[1])           ## define the dimension of data axes
y = np.linspace(0, 1, z.shape[0])  

## Load permittivity
h5file = h5py.File('structure.h5', "r")
data = np.array(h5file['eps'])[0,:,:]
eps = data

## Start figure + subplot (interactive)
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, axisbg='w', aspect=(float(len(y))/len(x)))
fig.subplots_adjust(left=.05, bottom=.05, right=.99, top=.99, wspace=.05, hspace=.05)

## Plot contours for gridded data
#extent = max(-np.min(z), np.max(z))
#contours = plt.contourf(x, y, z, levels=np.linspace(-extent,extent,50), cmap=matplotlib.cm.RdBu, extend='both')

if (len(sys.argv)>1) and ('angle'  in sys.argv[1]):
    contours = plt.contourf(x, y, np.angle(z), levels=np.linspace(-np.pi,np.pi,50), cmap=matplotlib.cm.hsv, extend='both')
else:
    contours = plt.contourf(x, y, np.log10(np.abs(z)), levels=np.linspace(-3,0,50), cmap=matplotlib.cm.gist_earth, extend='both')

for contour in contours.collections: contour.set_antialiased(False)     ## optional: avoid white aliasing (for matplotlib 1.0.1 and older) 
plt.contour(x, y, eps, levels=[3,3], colors='k',lw=.5)      ## optional: plot black contour at zero
#plt.colorbar()                                                          ## optional: colorbar


## ==== Outputting ====
## Finish the plot + save 
plt.xlabel(u"x"); 
plt.ylabel(u"y"); 
plt.grid()
plt.legend(prop={'size':10}, loc='upper right')
plt.savefig("output_%s_%s.png" % (dataset, '' if (len(sys.argv)==1) else sys.argv[1]), bbox_inches='tight')
