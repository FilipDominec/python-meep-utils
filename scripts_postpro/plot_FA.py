#!/usr/bin/env python
#-*- coding: utf-8 -*-

## ====== user settings =========
#timedomain = True
timedomain = False
c = 3e8
maxfreq = 2e12

## Import common moduli
import numpy as np
from scipy.constants import c, hbar, pi
import matplotlib, sys, os, time
import matplotlib.pyplot as plt

## Start figure + subplot (interactive)
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, axisbg='w')
fig.subplots_adjust(left=.05, bottom=.05, right=.99, top=.99, wspace=.05, hspace=.05)

## Decide the filename to load data
import sys, os
if len(sys.argv)>0: 
    params = dict([x.split('=')    for x in os.getcwd().split('_')   if '=' in x])
    simtime = float(params.get('simtime', 100e-12))
    size_y = float(params.get('size_y', 2200e-6))
    print 'simtime, size_y', simtime, size_y


filename = sys.argv[1] if len(sys.argv)>1 else 'SliceByVolume.h5'
if not os.path.isfile(filename): raise IOError, 'File %s can not be opened!' % filename

## Load n-dimensional arrays from a HDF5 file
import h5py
h5file = h5py.File(filename, "r")
print "Found datasets:", h5file.keys()
time1 = time.time()
data  = np.array(h5file['ex.r'])
data += np.array(h5file['ex.i']) * 1j
print "Loaded dataset with shape:", data.shape, 'in %04d s.' % (time.time()-time1)
try:
    Et = data[:,-1,:]       ## take the farthest slice by the z-axis
except IndexError:
    Et = data               ## if data already 2D
print simtime, Et.shape[1]
#print np.linspace(0, float(simtime), Et.shape[1])
t = np.linspace(0, simtime, Et.shape[1])           ## define the dimension of data axes
y = np.linspace(0, size_y, Et.shape[0])  

## Fourier transform along the time axis
freq    = np.fft.fftfreq(len(t), d=(t[1]-t[0]))         # calculate the frequency axis with proper spacing
Efy     = np.fft.fft(Et, axis=1) / len(t) * 2*np.pi     # calculate the FFT values

freq    = np.fft.fftshift(freq) #+ freq[len(freq)/2]
Efy     = np.fft.fftshift(Efy)

## Fourier transform along the spatial axis
kT      = np.fft.fftfreq(len(y), d=(y[1]-y[0]))        # calculate the frequency axis with proper spacing
Ef      = np.fft.fft(Efy, axis=0) / len(y) * 2*np.pi     # calculate the FFT values

kT      = np.fft.fftshift(kT)
def fftshift2(arr): return np.vstack([arr[len(arr)/2:,:], arr[:len(arr)/2,:]])
Ef =  fftshift2(Ef)

## Truncate too high frequencies
truncated = np.logical_and(freq>0, freq<maxfreq)         # (optional) get the frequency range
freq = freq[truncated]
Ef   = Ef[:,truncated]

## plot contours for gridded data
if timedomain: 
    #contours = plt.contourf(t, y, np.log10(np.abs(et)+1e-6), cmap=matplotlib.cm.gist_earth, extend='both') # levels=np.arange(0.,1,.01), 
    contours = plt.contourf(t, y, Et, cmap=matplotlib.cm.RdBu, extend='both') # levels=np.arange(0.,1,.01), 
else:
    toplot = np.log10(np.abs(Ef))
    contours = plt.contourf(freq, kT, toplot,  cmap=matplotlib.cm.gist_earth, levels=np.linspace(np.min(toplot)*.8+np.max(toplot)*.2,np.max(toplot),200) ,extend='both') #  
    #contours = plt.contourf(freq, kT, np.abs(Ef), cmap=matplotlib.cm.gist_earth, extend='both') # levels=np.arange(0.,1,.01), 
    plt.plot([0, maxfreq], [0, 0], c='w',lw=.5)
    plt.plot([0, maxfreq], [0, maxfreq/c], c='w',lw=.5)
    plt.plot([0, maxfreq], [0, -maxfreq/c], c='w',lw=.5)
    plt.annotate('+45$^\\circ$', xy = (maxfreq/2, maxfreq/2/c), xytext = (-10, 10), textcoords='offset points',color='w')
    plt.annotate('-45$^\\circ$', xy = (maxfreq/2, -maxfreq/2/c), xytext = (10, 10), textcoords='offset points',color='w')
for contour in contours.collections: contour.set_antialiased(False)     ## optional: avoid white aliasing (for matplotlib 1.0.1 and older) 
plt.colorbar()                                                          ## optional: colorbar

## Plot prior-calculated index of refraction
#try:
(f, N) = np.loadtxt(params.get('nref', './effparam.dat'), usecols=(0,5), unpack=True)
truncated = np.logical_and(f>0, f<maxfreq)         # (optional) get the frequency range
f = f[truncated]
N = N[truncated]

## Compute expected refraction
alpha = np.arctan(.5)      ## angle of the front face;  0.5 means the wedge front face is defined by "z - y*0.5 = const" 
chi   = alpha - np.arcsin(np.sin(alpha)/N)        ## internal angle in the wedge
delta = np.arcsin(N * np.sin(alpha - np.arcsin(np.sin(alpha)/N)))

## Plot lines for refraction
print N
plt.plot(f,  delta * f/3e8, color="#FF8800", label=u"$y'$", ls='-', c='k',lw=1)
plt.plot(f, -delta * f/3e8, color="#FF8800", label=u"$y'$", ls='-', c='k',lw=.2)

## Plot lines for N
plt.plot(f, np.real(N)*1000, color="#FF8800", label=u"$y'$", ls='-', c='w',lw=1)
plt.plot(f, np.imag(N)*1000, color="#FF8800", label=u"$y'$", ls='-', c='w',lw=1)
#except:
    #print "Warning: Refractive index reference could not be loaded/plot", sys.exc_info()[0]

## Finish the plot + save 
if not timedomain: plt.ylim((-2e4,2e4))
plt.xlabel(u"frequency"); 
plt.ylabel(u"tangential wavenumber K$_y$"); 
plt.grid()
plt.legend(prop={'size':10}, loc='upper right')
plt.savefig("output_%s.png" % ('_t' if timedomain else ''), bbox_inches='tight')


