#!/usr/bin/env python
#-*- coding: utf-8 -*-
## Import common moduli
import numpy as np
from scipy.constants import c, hbar, pi
import matplotlib, os, sys, argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--fieldfile',  type=str,                   help='HDF5 file containing the field record to be Fourier-transformed')
parser.add_argument('--nfile',      type=str,   default='',     help='optional data file containing the retrieved index of refraction for comparison')
parser.add_argument('--title',      type=str,   default='',     help='plot title')
parser.add_argument('--resolution', type=float, default=2e-6,   help='for realistic scaling of the spatial axis')
parser.add_argument('--cellnumber', type=int,   default=3,      help='number of cells') ## todo update code
parser.add_argument('--cellsize',   type=float, default=100e-6, help='z-dimensions of the simulated cells')
parser.add_argument('--simtime',    type=float, default=500e-12,help='time of simulations provided')
parser.add_argument('--frequnit',   type=float, default=1e12,   help='unit magnitude on the frequency axis')

parser.add_argument('--freqlim1',   type=str,   default='',     help='start for the plotted parameter range')
parser.add_argument('--freqlim2',   type=str,   default='',     help='end for the plotted parameter range')
parser.add_argument('--Nlim1',      type=str,   default='',     help='lower value of the effective index, i.e. both N.real N.imag')
parser.add_argument('--Nlim2',      type=str,   default='',     help='upper value of the effective index, i.e. both N.real N.imag')
parser.add_argument('--fieldlim1',  type=str,   default='',     help='lower contour value for the field')
parser.add_argument('--fieldlim2',  type=str,   default='',     help='upper contour value for the field')
parser.add_argument('--fieldlabel', type=str,   default='Field amplitude $|\mathbf E|$', help='just a text over the field contours')
parser.add_argument('--Nconj',      type=str,   default='',     help='set to "yes" to plot complex conjugated values') 
parser.add_argument('--logfield',   type=str,   default='yes',  help='plot the field as logarithmic')
parser.add_argument('--colormap',   type=str,   default='default', help='matplotlib colormap, available are: gist_earth (default for contours), jet, hsv, greys, dark2, brg...')
parser.add_argument('--numcontours',type=int,   default=30,     help='number of levels in the contour plot (default 50)')
parser.add_argument('--contourresx',type=int,   default=200,    help='row length of the internal interpolation matrix for contour plot (default 200)')
parser.add_argument('--contourresp',type=int,   default=200,    help='column height of the internal interpolation matrix for contour plot (default 200)')
parser.add_argument('--usetex',     type=str,   default='yes', help='by default, LaTeX is used for nicer typesetting')
args = parser.parse_args()

print args.frequnit
print args.freqlim1
print args.freqlim2

## Plotting style
if args.usetex.lower() in ('yes', 'true'): 
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{txfonts}, \usepackage{upgreek}') #, \usepackage{palatino} \usepackage{lmodern}, 
matplotlib.rc('font', size=12)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})  ## select fonts

cmap = matplotlib.cm.gist_earth if (args.colormap == 'default') else getattr(matplotlib.cm, args.colormap)  

##Start figure + subplot
fig = plt.figure(figsize=(6,7))   # for publication
ax = plt.subplot(211, axisbg='w')
fig.subplots_adjust(left=.10, bottom=.08, right=.99, top=.99, wspace=.05, hspace=.05)

import h5py
try: 
    h5file = h5py.File(args.fieldfile, "r");  print 'File contains datasets:', h5file.keys()
    data  =  np.array((h5file[h5file.keys()[0]]))
except:
    print "Could not open file %s" % args.fieldfile; quit()

## select one field stripe to be processed (may average few of them)
Et = np.mean(data[0:1],axis=0)
t = np.linspace(0, args.simtime, Et.shape[1])

## convert time-axis to frequency-axis
freq = np.fft.fftshift(np.fft.fftfreq(len(t), d=(t[1]-t[0])))       # calculate the frequency axis with proper spacing
Ef = np.zeros_like(Et) * (1+0j)                                     # prepare a complex array
for m in range(Ef.shape[0]):                                        ## go through all z-positions
    Ef[m] = np.fft.fftshift(np.fft.fft(Et[m])) / len(Et[0]) * 2*np.pi                # calculate the FFT values, normalize to number of data points

## Normalize the pseudo-gaussian source amplitude
#Ef = Ef / (np.exp(-((freq-2.5e12)/.8e12)**2) + np.exp(-((freq+2.5e12)/.8e12)**2))

## Get the simulation dimension along the z-axis
z_range = args.resolution * data.shape[1] / args.cellsize
z_cellofs = (z_range - args.cellnumber)/2
#print 'z_range=', z_range, '*args.cellsize; ', 'z_cellofs', z_cellofs
z_axis = np.linspace(0, z_range, Ef.shape[0])
z_axis -= z_cellofs

## Plot contours for gridded data
field_to_plot = np.log10(np.abs(Ef)) if (args.logfield.lower() == 'yes') else Ef 
fieldlim1 = np.min(Ef) if args.fieldlim1 == "" else float(args.fieldlim1)
fieldlim2 = np.max(Ef) if args.fieldlim2 == "" else float(args.fieldlim2)
print " fieldlim1, fieldlim2", fieldlim1, fieldlim2
CS = plt.contourf(freq/args.frequnit, z_axis, field_to_plot,
        levels=np.linspace(fieldlim1, fieldlim2, args.numcontours), cmap=cmap, extend='both')
for contour in CS.collections: contour.set_antialiased(False)
if args.fieldlabel != "": plt.plot([],[], lw=0, label=args.fieldlabel)

#plt.colorbar()
#plt.contour(freq, z_axis, Ef.real, levels=[0,0], colors='k', extend='both') ## Shows the phase
#plt.contour(freq, z_axis, Ef.imag, levels=[0,0], colors='w', extend='both')
#plt.contour(freq, z_axis, np.angle(Ef), levels=[0,0], colors='k', extend='both') ## Shows the phase

## Plot complex amplitude
#plt.contourf(freq, z_axis, np.angle(Ef+1j*Ef[::-1,:]), levels=np.linspace(-np.pi, np.pi), cmap=matplotlib.cm.hsv, extend='both')
#plt.contourf(freq, z_axis, np.log(np.abs(Ef)), levels=np.linspace(-15, np.max(np.log(np.abs(Ef)))), cmap=matplotlib.cm.gray, extend='both', alpha=.9)

## Plot cell boundaries
for zpos in range(args.cellnumber+1):
    plt.plot([0, np.max(freq)], [zpos,zpos], c='w', lw=1, scaley=False)

plt.ylabel('Position on the $z$-axis [$a$]')

bbox        = dict(boxstyle='round, pad=.15', fc='white', alpha=1)
arrowprops  = dict(arrowstyle = '<->', color='w')
plt.annotate('', xy = (420e9/args.frequnit, 4), xytext = (1020e9/args.frequnit, 4),
        textcoords='data', ha='center', va='bottom', bbox=bbox, arrowprops=arrowprops)
plt.annotate('band gap', xy = (720e9/args.frequnit, 4), xytext = (720e9/args.frequnit, 4),
        textcoords='data', ha='center', va='bottom', bbox=None, arrowprops=None, color='w')

plt.annotate('', xy = (1050e9/args.frequnit, 4), xytext = (1150e9/args.frequnit, 4),
        textcoords='data', ha='center', va='bottom', bbox=bbox, arrowprops=arrowprops)
plt.annotate('band gap', xy = (1050e9/args.frequnit, 4), xytext = (1100e9/args.frequnit, 4),
        textcoords='data', ha='center', va='bottom', bbox=None, arrowprops=None, color='w')

plt.annotate('1$^{st}$ cell', xy = (220e9/args.frequnit, .5), xytext = (200e9/args.frequnit, .5),
        textcoords='data', ha='center', va='center', bbox=bbox, arrowprops=None)
plt.annotate('2$^{nd}$ cell', xy = (220e9/args.frequnit, 1.5), xytext = (200e9/args.frequnit, 1.5),
        textcoords='data', ha='center', va='center', bbox=bbox, arrowprops=None)
plt.annotate('3$^{rd}$ cell', xy = (220e9/args.frequnit, 2.5), xytext = (200e9/args.frequnit, 2.5),
        textcoords='data', ha='center', va='center', bbox=bbox, arrowprops=None)

plt.annotate('', xy = (800e9/args.frequnit, 1.3), xytext = (20,20),
        textcoords='offset points', ha='center', va='bottom', bbox=bbox, arrowprops=dict(arrowstyle = 'simple', color='k'))
plt.annotate('', xy = (950e9/args.frequnit, 1.3), xytext = (20,20),
        textcoords='offset points', ha='center', va='bottom', bbox=bbox, arrowprops=dict(arrowstyle = 'simple', color='k'))

plt.legend(loc='lower left')


if args.nfile != '':
    ax.label_outer()
    ax = plt.subplot(212, axisbg='w', sharex=ax)
    (f, Nr, Ni) = np.loadtxt(args.nfile, usecols=[0,5,6], unpack=True)
    Kr = Nr*f*np.pi*2*args.cellsize/c
    Ki = Ni*f*np.pi*2*args.cellsize/c
    zero_axis = np.zeros_like(f)+0
    n_scaling = 1
    plt.plot(f/args.frequnit,(Nr+zero_axis)/n_scaling, c='b', lw=2, label="real part $N'$")
    plt.plot(f/args.frequnit,(Ni+zero_axis)/n_scaling, c='b', lw=2,ls='--', label="imaginary part $N''$")
    plt.legend(loc='lower left')
    #plt.plot(f/args.frequnit,(Kr+zero_axis)/n_scaling, c='k', lw=2)
    #plt.plot(f/args.frequnit,(Ki+zero_axis)/n_scaling, c='k', lw=2,ls='--')
    #plt.plot(f,zero_axis/n_scaling, c='k', lw=2,ls=':')
    #plt.plot(f,(zero_axis+np.ones_like(f))/n_scaling, c='k', lw=2, ls=':')

    #plt.annotate('', xy = (800e9/args.frequnit, 1.0), xytext = (20,20),
            #textcoords='offset points', ha='center', va='bottom', bbox=bbox, arrowprops=dict(arrowstyle = 'simple', color='k'))
    #plt.annotate('', xy = (950e9/args.frequnit, -0.8), xytext = (20,20),
            #textcoords='offset points', ha='center', va='bottom', bbox=bbox, arrowprops=dict(arrowstyle = 'simple', color='k'))

    plt.grid()
    plt.ylabel('Refractive index $N_{\\text{eff}}$')

if args.freqlim1 != '': plt.xlim(left=float(args.freqlim1)/args.frequnit)
if args.freqlim2 != '': plt.xlim(right=float(args.freqlim2)/args.frequnit) 
# TODO ylim
plt.xlabel('Frequency [THz]')
plt.title(args.title)
    
plt.savefig("fd_ampli_publix.png")
