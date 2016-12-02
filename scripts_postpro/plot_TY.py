#!/usr/bin/env python
#-*- coding: utf-8 -*-

simtime = 80e-12
size_y = 1400e-6
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


## Start figure + subplot (interactive)
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111, axisbg='w')
fig.subplots_adjust(left=.05, bottom=.05, right=.99, top=.99, wspace=.05, hspace=.05)

## Decide the filename to load data
import sys
filename = sys.argv[1] if len(sys.argv)>1 else 'input.dat'
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
t = np.linspace(0, simtime, Et.shape[1])           ## define the dimension of data axes
y = np.linspace(0, size_y, Et.shape[0])  

## Export n-dimensional arrays to a HDF5 file


## Fourier transform
freq    = np.fft.fftfreq(len(t), d=(t[1]-t[0]))         # calculate the frequency axis with proper spacing
Efy     = np.fft.fft(Et, axis=1) / len(t) * 2*np.pi     # calculate the FFT values
#def ffts(arr):
    #return np.hstack([arr[len(arr)/2+1:], arr[:len(arr)/2]])
def ffts2(arr):
    return np.vstack([arr[len(arr)/2:,:], arr[:len(arr)/2,:]])
#freq = ffts(freq)
#Efy =  ffts2(Efy)

freq    = np.fft.fftshift(freq) #+ freq[len(freq)/2]
Efy     = np.fft.fftshift(Efy)

kT      = np.fft.fftfreq(len(y), d=(y[1]-y[0]))        # calculate the frequency axis with proper spacing
Ef      = np.fft.fft(Efy, axis=0) / len(y) * 2*np.pi     # calculate the FFT values

kT      = np.fft.fftshift(kT)
#Ef      = np.fft.fftshift(Ef)
print Ef.shape
Ef =  ffts2(Ef)
print Ef.shape

truncated = np.logical_and(freq>0, freq<maxfreq)         # (optional) get the frequency range
freq = freq[truncated]
Ef   = Ef[:,truncated]

print 'freq', freq.shape, freq[::10]
print 'kT', kT.shape, kT[::10]


## plot contours for gridded data
#contours = plt.contourf(t, y, np.log10(np.abs(et)+1e-6), cmap=matplotlib.cm.gist_earth, extend='both') # levels=np.arange(0.,1,.01), 
#contours = plt.contourf(t, y, et, cmap=matplotlib.cm.rdbu, extend='both') # levels=np.arange(0.,1,.01), 

toplot = (np.abs(Et))
contours = plt.contourf(t, y, toplot,  cmap=matplotlib.cm.gist_earth, levels=np.linspace(np.min(toplot)*0+np.max(toplot)*0,np.max(toplot),200) ,extend='both') #  
#contours = plt.contourf(freq, kT, np.abs(Ef), cmap=matplotlib.cm.gist_earth, extend='both') # levels=np.arange(0.,1,.01), 

#plt.plot([0, maxfreq], [0, 0], c='w',lw=.5)
#plt.plot([0, maxfreq], [0, maxfreq/c], c='w',lw=.5)
#plt.plot([0, maxfreq], [0, -maxfreq/c], c='w',lw=.5)
#plt.annotate('+45$^\\circ$', xy = (maxfreq/2, maxfreq/2/c), xytext = (-10, 10), textcoords='offset points',color='w')
#plt.annotate('-45$^\\circ$', xy = (maxfreq/2, -maxfreq/2/c), xytext = (10, 10), textcoords='offset points',color='w')
#

try:
    ## Load 1D curve
    filename = "effparam.dat"
    (x, y) = np.loadtxt(filename, usecols=(0,5), unpack=True)
    truncated = np.logical_and(x>0, x<maxfreq)         # (optional) get the frequency range
    x = x[truncated]
    y = y[truncated]
    ## Plot line
    plt.plot(x, np.real(y)*1000, color="#FF8800", label=u"$y'$", ls='-', c='w',lw=1)

except:
    print "refractive index could not be loaded"


for contour in contours.collections: contour.set_antialiased(False)     ## optional: avoid white aliasing (for matplotlib 1.0.1 and older) 
plt.colorbar()                                                          ## optional: colorbar

## Finish the plot + save 
#plt.ylim((-2e4,2e4))
plt.xlabel(u"time"); 
plt.ylabel(u"y"); 
plt.grid()
plt.legend(prop={'size':10}, loc='upper right')
plt.savefig("output_T-Y.png", bbox_inches='tight')


