#!/usr/bin/env python
#-*- coding: utf-8 -*-
## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, hbar, pi

scale = 250  ## scale of vector glyphs
decimate = 1
zoom = 2.5  ## ratio of full simulation width to the shown area
resolution = 3e-6
sizeunit = 1e-6
xcenter = resolution/3
ycenter = resolution/2 *0

## Use LaTeX
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', size=12)
matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{yfonts}, \usepackage{txfonts}, \usepackage{lmodern},')
## Start figure + subplot (interactive)
fig = plt.figure(figsize=(4,3))
ax = plt.subplot(111, axisbg='w', adjustable='box', aspect=1.0)
fig.subplots_adjust(left=.05, bottom=.05, right=.99, top=.99, wspace=.05, hspace=.05)

## Load from file
def load_hdf(filename, dataset):
    import h5py
    h5file = h5py.File(filename, "r")
    data = np.array(h5file[dataset])
    ysize, zsize = data.shape[1], data.shape[2]
    fromz, toz = int((zsize-ysize/zoom)/2)+1, int((zsize+ysize/zoom)/2) 
    fromy, toy = int((ysize-ysize/zoom)/2)+1,   int((ysize+ysize/zoom)/2) 
    data = data[int(data.shape[0]/2+0),fromy:toy,fromz:toz,0]
    return data

Ex = load_hdf('snapshot_Ex_t0.000000e+00.h5', 'ex.r')
Hy = load_hdf('snapshot_Hy_t0.000000e+00.h5', 'hy.r')
Hz = load_hdf('snapshot_Hz_t0.000000e+00.h5', 'hz.r')
## Generate axes
#y = np.arange(-Ex.shape[1]*resolution/sizeunit/2, Ex.shape[1]*resolution/sizeunit/2, resolution/sizeunit)           ## define the dimension of data axes
#z = np.arange(-Ex.shape[0]*resolution/sizeunit/2, Ex.shape[0]*resolution/sizeunit/2, resolution/sizeunit)           ## define the dimension of data axes
y = np.linspace(-.5, .5, Ex.shape[1]) * Ex.shape[1]*resolution/sizeunit - xcenter/sizeunit         ## define the dimension of data axes
z = np.linspace(-.5, .5, Ex.shape[0]) * Ex.shape[0]*resolution/sizeunit - ycenter/sizeunit          ## define the dimension of data axes
print len(y), len(z)
ys, zs = np.meshgrid(y,z)

## Plot contours for gridded data
#extent = max(-np.min(Ex), np.max(Ex))  ## for balanced color scale 
extent = 1
Ex = Ex / 1.2
contours = plt.contourf(y, z, Ex, levels=np.linspace(-extent,extent,50), cmap=matplotlib.cm.RdBu_r, extend='both')
for contour in contours.collections: contour.set_antialiased(False)     ## optional: avoid white aliasing (for matplotlib 1.0.1 and older) 
plt.colorbar().set_ticks([-1, 0, 1])                                                     ## optional: colorbar

#ys, zs, Hz, Hy = [arr[::decimate,::decimate] for arr in (ys, zs, Hz, Hy)]
#q = plt.quiver(ys, zs, Hz, Hy, linewidths=(0,), edgecolors=('k'), 
        #headwidth=7, headlength=11, headaxislength=9, pivot='mid', scale=scale)
#plt.quiverkey(q, .1, .1, 2, 'magnetic field')

angle = np.linspace(0, 2*np.pi,100)
plt.plot((np.sin(angle)*12.5e-6)/sizeunit, (np.cos(angle)*12.5e-6)/sizeunit, c='k', ls='--')


plt.xlabel(u"z [$\muup$m]"); 
plt.ylabel(u"y [$\muup$m]"); 
plt.grid()
plt.legend(prop={'size':10}, loc='upper right')
plt.savefig("output_new.pdf", bbox_inches='tight')
