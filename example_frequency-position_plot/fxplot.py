#!/usr/bin/env python
#-*- coding: utf-8 -*-
## Import common moduli
import numpy as np
from scipy.constants import c, hbar, pi
import matplotlib, os, sys
import matplotlib.pyplot as plt
matplotlib.rc('text', usetex=True)
matplotlib.rc('font', size=14)
#matplotlib.rc('text.latex', preamble = \
        #'\usepackage{amsmath}, ')

matplotlib.rc('font',**{'family':'sans-serif', 'sans-serif':['Computer Modern Sans serif']})
##Start figure + subplot (interactive)
#fig = plt.figure(figsize=(19.2,10.8))   # for full-screen
fig = plt.figure(figsize=(6,7))   # for publication
ax = plt.subplot(211, axisbg='w')
fig.subplots_adjust(left=.10, bottom=.08, right=.99, top=.99, wspace=.05, hspace=.05)

resolution = 2e-6
cells_n = 3
cellspacing = 100e-6
t_range = 500e-12
frequnit = 1e12
print "t_range=", t_range


import h5py
h5file = h5py.File('Slice_x0.0000.h5', "r");  print 'File contains datasets:', h5file.keys()
data  =   np.array((h5file[h5file.keys()[0]]))

## select one y-stripe to be processed (may average few of them, not important?)
Et = np.mean(data[0:1],axis=0)
t = np.linspace(0, t_range, Et.shape[1])

## convert time-axis to frequency-axis
freq = np.fft.fftshift(np.fft.fftfreq(len(t), d=(t[1]-t[0])))       # calculate the frequency axis with proper spacing
Ef = np.zeros_like(Et) * (1+0j)                                     # prepare a complex array
for n in range(Ef.shape[0]):                                        ## go through all z-positions
    Ef[n] = np.fft.fftshift(np.fft.fft(Et[n])) / len(Et[0]) * 2*np.pi                # calculate the FFT values, normalize to number of data points

## Normalize the pseudo-gaussian source amplitude
Ef = Ef / (np.exp(-((freq-2.5e12)/.8e12)**2) + np.exp(-((freq+2.5e12)/.8e12)**2))

## Get the simulation dimension along the z-axis
z_range = resolution * data.shape[1] / cellspacing
z_cellofs = (z_range - cells_n)/2
#print 'z_range=', z_range, '*cellspacing;; ', 'z_cellofs', z_cellofs
z_axis = np.linspace(0, z_range, Ef.shape[0])
z_axis -= z_cellofs

## Plot contours for gridded data
CS = plt.contourf(freq/frequnit, z_axis, np.log(np.abs(Ef)), levels=np.linspace(-15,-5), cmap=matplotlib.cm.gist_earth, extend='both')
for contour in CS.collections: contour.set_antialiased(False)
plt.plot([],[], lw=0, label='Field amplitude $|\mathbf E|$')

#plt.colorbar()
#plt.contour(freq, z_axis, Ef.real, levels=[0,0], colors='k', extend='both') ## Shows the phase
#plt.contour(freq, z_axis, Ef.imag, levels=[0,0], colors='w', extend='both')
#plt.contour(freq, z_axis, np.angle(Ef), levels=[0,0], colors='k', extend='both') ## Shows the phase

## Plot complex amplitude
#plt.contourf(freq, z_axis, np.angle(Ef+1j*Ef[::-1,:]), levels=np.linspace(-np.pi, np.pi), cmap=matplotlib.cm.hsv, extend='both')
#plt.contourf(freq, z_axis, np.log(np.abs(Ef)), levels=np.linspace(-15, np.max(np.log(np.abs(Ef)))), cmap=matplotlib.cm.gray, extend='both', alpha=.9)

## Plot cell boundaries
for zpos in range(cells_n+1):
    plt.plot([0, np.max(freq)], [zpos,zpos], c='w', lw=1)

plt.ylabel('Position on the $z$-axis [$a$]')

bbox        = dict(boxstyle='round, pad=.15', fc='white', alpha=1)
arrowprops  = dict(arrowstyle = '<->', color='w')
plt.annotate('', xy = (420e9/frequnit, 4), xytext = (1020e9/frequnit, 4),
        textcoords='data', ha='center', va='bottom', bbox=bbox, arrowprops=arrowprops)
plt.annotate('band gap', xy = (720e9/frequnit, 4), xytext = (720e9/frequnit, 4),
        textcoords='data', ha='center', va='bottom', bbox=None, arrowprops=None, color='w')

plt.annotate('', xy = (1050e9/frequnit, 4), xytext = (1150e9/frequnit, 4),
        textcoords='data', ha='center', va='bottom', bbox=bbox, arrowprops=arrowprops)
plt.annotate('band gap', xy = (1050e9/frequnit, 4), xytext = (1100e9/frequnit, 4),
        textcoords='data', ha='center', va='bottom', bbox=None, arrowprops=None, color='w')

plt.annotate('1$^{st}$ cell', xy = (220e9/frequnit, .5), xytext = (200e9/frequnit, .5),
        textcoords='data', ha='center', va='center', bbox=bbox, arrowprops=None)
plt.annotate('2$^{nd}$ cell', xy = (220e9/frequnit, 1.5), xytext = (200e9/frequnit, 1.5),
        textcoords='data', ha='center', va='center', bbox=bbox, arrowprops=None)
plt.annotate('3$^{rd}$ cell', xy = (220e9/frequnit, 2.5), xytext = (200e9/frequnit, 2.5),
        textcoords='data', ha='center', va='center', bbox=bbox, arrowprops=None)

plt.annotate('', xy = (800e9/frequnit, 1.3), xytext = (20,20),
        textcoords='offset points', ha='center', va='bottom', bbox=bbox, arrowprops=dict(arrowstyle = 'simple', color='k'))
plt.annotate('', xy = (950e9/frequnit, 1.3), xytext = (20,20),
        textcoords='offset points', ha='center', va='bottom', bbox=bbox, arrowprops=dict(arrowstyle = 'simple', color='k'))

plt.legend(loc='lower left')


if len(sys.argv) > 1:
    ax.label_outer()
    ax = plt.subplot(212, axisbg='w', sharex=ax)
    (f, Nr, Ni) = np.loadtxt(sys.argv[1], usecols=[0,5,6], unpack=True)
    Kr = Nr*f*np.pi*2*cellspacing/3e8
    Ki = Ni*f*np.pi*2*cellspacing/3e8
    zero_axis = np.zeros_like(f)+0
    n_scaling = 1
    plt.plot(f/frequnit,(Nr+zero_axis)/n_scaling, c='b', lw=2, label="real part $N'$")
    plt.plot(f/frequnit,(Ni+zero_axis)/n_scaling, c='b', lw=2,ls='--', label="imaginary part $N''$")
    plt.legend(loc='lower left')
    #plt.plot(f/frequnit,(Kr+zero_axis)/n_scaling, c='k', lw=2)
    #plt.plot(f/frequnit,(Ki+zero_axis)/n_scaling, c='k', lw=2,ls='--')
    #plt.plot(f,zero_axis/n_scaling, c='k', lw=2,ls=':')
    #plt.plot(f,(zero_axis+np.ones_like(f))/n_scaling, c='k', lw=2, ls=':')

    plt.annotate('', xy = (800e9/frequnit, 1.0), xytext = (20,20),
            textcoords='offset points', ha='center', va='bottom', bbox=bbox, arrowprops=dict(arrowstyle = 'simple', color='k'))
    plt.annotate('', xy = (950e9/frequnit, -0.8), xytext = (20,20),
            textcoords='offset points', ha='center', va='bottom', bbox=bbox, arrowprops=dict(arrowstyle = 'simple', color='k'))

    plt.grid()
    plt.ylabel('Refractive index $N$')

#plt.ylim((0,z_range))
#plt.ylim((0.2*z_range,0.6*z_range))
plt.xlim((200e9/frequnit, 1400e9/frequnit))
plt.xlabel('Frequency [THz]')
#plt.xlim((0, min(np.max(freq), 1.75e12)))
    
    #plt.show()
#plt.savefig("tdc.png", bbox_inches='tight')
plt.savefig("fd_ampli_publix.png")
#plt.savefig("fd_phase.png", bbox_inches='tight')
#plt.savefig("fd_complex.png", bbox_inches='tight')
