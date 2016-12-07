#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Import common moduli
import matplotlib, sys, os, time, argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, hbar, pi

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--xlim1',      type=str,   default='', help='start for the x-axis range')
parser.add_argument('--xlim2',      type=str,   default='', help='end for the x-axis range')
parser.add_argument('--ylim1',      type=str,   default='0.', help='start for the plotted value range')
parser.add_argument('--ylim2',      type=str,   default='', help='end for the plotteld value range')
parser.add_argument('--numcontours',type=int,   default=50, help='number of levels in the contour plot (default 50)')
parser.add_argument('--contourresx',type=int,   default=200,help='row length of the internal interpolation matrix for contour plot (default 200)')
parser.add_argument('--contourresp',type=int,   default=200,help='column height of the internal interpolation matrix for contour plot (default 200)')
parser.add_argument('filenames',    type=str,   nargs='+', help='CSV files to be processed')
#parser.add_argument('--output',     type=str,   default='*.png', help='output file; *.png or *.pdf (etc.) auto-selects a name with the given format')
args = parser.parse_args()

## Use LaTeX
matplotlib.rc('text', usetex=True)
matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{txfonts}, \usepackage{upgreek}') 
# optionally: \usepackage{palatino} \usepackage{lmodern}, 
matplotlib.rc('font', size=12)
matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})  ## select fonts

# -- settings --
minfreq = .5e12
maxfreq = 1e12 
frequnit = 1

plot_FFT = True
use_vacuum_ref = False
interp_anisotropy = 2e-5    # value lower than 1. interpolates rather vertically; optimize if plot disintegrates
FFTcutoff = 0.8             # Hann-like window to suppress spectral leakage in FFT (mostly for aesthetic purposes)

plot_FDM = True
FDMtrunc = (.1, .5)         # Harminv (also known as FDM) may work better when it is supplied a shorter time record
                            # and it requires clipping the initial timespan when the source is operating

plot_NRef = True


cmap = matplotlib.cm.Blues_r
cmap = matplotlib.colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=.4, b=1.0), cmap(np.linspace(.4, 1.0, 100)))


if   '-phase' in sys.argv[1:]: 
    filesuffix = 'phase'
elif '-epsilon' in sys.argv[1:]: 
    filesuffix = 'epsilon'
else:
    filesuffix = 'ampli'


## Load and prepare data (from multiple files)
#filenames = args.filenames
#if len(filenames) == 0: print "Error: no data file to be plotted was provided as argument" ; quit()
Efs = []
EfsRef = []
Kzs = []     ## TODO allow also scanning over Kx, Ky (which allows for plotting dispersion curves also along "Γ-M" and other directions in K-space)
freqs = []
fs = []

FDM_freqs = []
FDM_amplis = []
FDM_phases = []
FDM_Kzs = []

for vacuum_ref in [True, False] if use_vacuum_ref else [False]:
    if vacuum_ref:  filenames = os.listdir('ref')
    else:           filenames = args.filenames
    for filename, color in zip(filenames, matplotlib.cm.hsv(np.linspace(0,1,len(filenames)))): 
        Kz = None
        cellsize = None
        with open(filename) as datafile:
            for line in datafile:
                if line[0:1] in "0123456789": break         # end of file header
                value = line.replace(",", " ").split()[-1]  # the value of the parameter will be separated by space or comma
                if not Kz and ("Kz" in line): Kz = float(value)
                if not cellsize and ("cellsize" in line): cellsize = float(value)
        Kz = Kz*cellsize/2/np.pi
        (t, E) = np.loadtxt(filename, usecols=list(range(2)), unpack=True, )

        if plot_FDM and not vacuum_ref:
            import harminv_wrapper
            tscale = 3e9            ## TODO check again that this prescaling is needed
            t1 = t[len(t)*FDMtrunc[0]:len(t)*FDMtrunc[1]]*tscale
            t1 -= np.min(t1)
            E1 = E[len(t)*FDMtrunc[0]:len(t)*FDMtrunc[1]]
            try:
                hi = harminv_wrapper.harminv(t1, E1, d=.1, f=15)
                hi['frequency'] *= tscale /frequnit  * cellsize/3e8/2/np.pi
                hi['amplitude'] /= np.max(hi['amplitude'])
                hi['error'] /= np.max(hi['amplitude'])
                FDM_freqs = np.append(FDM_freqs,  hi['frequency'])
                FDM_amplis= np.append(FDM_amplis,  hi['amplitude'])
                FDM_phases= np.append(FDM_phases,  hi['phase'])
                FDM_Kzs   = np.append(FDM_Kzs,    Kz*np.ones_like(hi['frequency'])) 
            except: 
                print "Error: Harminv did not find any oscillator in %s" % filename


        if plot_FFT:
            for field in (E,):
                field[t>max(t)*FFTcutoff] = field[t>max(t)*FFTcutoff]*(.5 + .5*np.cos(np.pi * (t[t>max(t)*FFTcutoff]/max(t)-FFTcutoff)/(1-FFTcutoff)))
            ## 1D FFT with cropping for useful frequencies
            freq    = np.fft.fftfreq(len(t), d=(t[1]-t[0]))  * cellsize/c/2/np.pi       # calculate the frequency axis with proper spacing
            Ef      = np.fft.fft(E, axis=0) / len(t) * 2*np.pi     # calculate the FFT values
            truncated = np.logical_and(freq>(minfreq*cellsize/3e8/2/np.pi), freq<(maxfreq*cellsize/3e8/2/np.pi))         # (optional) get the frequency range
            (Ef, freq) = map(lambda x: x[truncated], (Ef, freq))    # (optional) truncate the data points

            freq    = np.fft.fftshift(freq)
            pulsedelay  = 9.2e-12 ## obsolete
            Y      = np.fft.fftshift(Ef) / np.exp(-1j*pulsedelay*freq) * 100
            if not vacuum_ref:
                Ef = Y
                Efs     = np.append(Efs,    Ef)

                Kzs     = np.append(Kzs,    Kz*np.ones_like(freq))
                freqs   = np.append(freqs,  freq)
            else:
                EfsRef  = np.append(Efs,    Y)
                ## TODO: implement retrieval of the effective spatial-dispersive permittivity ε(ω,K)
                # if [sqrt(mu_r eps) c k / omega] - 1 = Y
                #Ef = ((freq*2*np.pi) / (Kz * c))**2 * (Y + 1)**2

                # if [sqrt(mu_r eps) c k / omega] - 1 = Y omega eps
                # then   a * q**2  +  b * q  +  c = 0
                # Where
                #omega = freq * 2*np.pi
                #A = 1/omega
                #B = c * Kz / omega**2
                #C = Y
                #q1 = -B  +  (B**2 - 4 * A * C)**.5
                #Ef = q1 ** (-1./2)
                #print Ef




## Plotting
plt.figure(figsize=(4,8))
ax = plt.subplot(111)
## Interpolate 2D grid from scattered data
from matplotlib.mlab import griddata
fi = np.linspace(0, maxfreq * cellsize/c/2/np.pi, 200)
ki = np.linspace(0, np.max(Kzs), 50)

## Plot contours for gridded data
if use_vacuum_ref:
    contour_data = Efs/EfsRef
else:
    contour_data = Efs

z = griddata(Kzs*interp_anisotropy, freqs/frequnit, np.abs(contour_data), ki*interp_anisotropy, fi, interp='linear')
log_min, log_max = np.log10(np.min(z)), np.log10(np.max(z))
log_min, log_max = log_min, log_min*.3 + log_max*.7
#log_min, log_max = -5, .5
levels = 10**(np.linspace(         log_min,          log_max,   args.numcontours))       ## where the contours are drawn
ticks  = 10**(np.arange(np.floor(log_min), np.ceil(log_max),  1))       ## where a number is printed
if plot_FFT:
    contours = plt.contourf(ki, fi, z, levels=levels, cmap=cmap, norm = matplotlib.colors.LogNorm())
    #plt.colorbar(ticks=ticks).set_ticklabels(['$10^{%d}$' % int(np.log10(t)) for t in ticks])

if plot_FFT:
    for contour in contours.collections: contour.set_antialiased(False)     ## optional: avoid white aliasing (for matplotlib 1.0.1 and older) 


if plot_NRef:
    try:
        f, Nre = np.loadtxt('NRef.dat', usecols=(0,5), unpack=True)
        plt.plot( Nre*f/c *cellsize, f * cellsize/c/2/np.pi, color='g', lw=1.5)               # positive-direction branches
        plt.plot(-Nre*f/c *cellsize, f * cellsize/c/2/np.pi, color='g', ls='--', lw=1.5)      # negative-direction branches
    except IOError:
        print "File NRef.dat was not found - not plotting the curve for comparison"


if plot_FDM:
    plt.scatter(FDM_Kzs, FDM_freqs, s=FDM_amplis*30+1, c=FDM_amplis) #, c=FDM_phases, cmap=plt.cm.hsv
 

## Simple axes
plt.ylim(ymin=minfreq * cellsize/c/2/np.pi, ymax=maxfreq * cellsize/c/2/np.pi) 
if args.ylim2 != "": plt.ylim(ymax=float(args.ylim2))
plt.xlim((np.min(ki), np.max(ki)))

plt.grid(True)


## Here we define the desired nonlinear dependence between dual y-axes:
def right_tick_function(p): return p/(cellsize/c/2/np.pi)/1e12
ax2 = ax.twinx()
#ax2.axis['right'].major_ticklabels.set_visible(False)
ax2.set_ylim(np.array(ax.get_ylim()))
ax2.set_ylabel('electron penetration depth (nm)')
## If we wish nice round numbers on the secondary axis ticks...
right_ax_limits = sorted([right_tick_function(lim) for lim in ax.get_ylim()])
yticks2 = matplotlib.ticker.MaxNLocator(nbins=8, steps=[1,2,5]).tick_values(*right_ax_limits)
## ... we must give them correct positions. To do so, we need to numerically invert the tick values:
from scipy.optimize import brentq
def right_tick_function_inv(r, target_value=0): 
    return brentq(lambda r,q: right_tick_function(r) - target_value, ax.get_ylim()[0], ax.get_ylim()[1], 1e6)
valid_ytick2loc, valid_ytick2val = [], []
for ytick2 in yticks2:
    try:
        valid_ytick2loc.append(right_tick_function_inv(0, ytick2))
        valid_ytick2val.append(ytick2)
    except ValueError:      ## (skip tick if the ticker.MaxNLocator gave invalid target_value to brentq optimization)
        pass
## Finally, we set the positions and tick values
ax2.set_yticks(valid_ytick2loc)
from matplotlib.ticker import FixedFormatter
ax2.yaxis.set_major_formatter(FixedFormatter(["%g" % righttick for righttick in valid_ytick2val]))


## The same approach for the x-axes (relation between photon energy and vavelength)
def top_tick_function(p): return p/(cellsize/2/np.pi) / 100
ax3 = ax.twiny()
#ax2.axis['top'].major_ticklabels.set_visible(False)
ax3.set_xlim(np.array(ax.get_xlim()))
ax3.set_xlabel('photon wavelength (nm)')
## If we wish nice round numbers on the secondary axis ticks...
print(ax.get_xlim())
print(top_tick_function(1.0))
print([top_tick_function(lim) for lim in ax.get_xlim()])
top_ax_limits = sorted([top_tick_function(lim) for lim in ax.get_xlim()])
print(top_ax_limits)
xticks2 = matplotlib.ticker.MaxNLocator(nbins=8, steps=[1,2,5]).tick_values(*top_ax_limits)
## ... we must give them correct positions. To do so, we need to numerically invert the tick values:
from scipy.optimize import brentq
def top_tick_function_inv(r, target_value=0): 
    return brentq(lambda r,q: top_tick_function(r) - target_value, ax.get_xlim()[0], ax.get_xlim()[1], 1e6)
valid_xtick2loc, valid_xtick2val = [], []
for xtick2 in xticks2:
    try:
        valid_xtick2loc.append(top_tick_function_inv(0, xtick2))
        valid_xtick2val.append(xtick2)
    except ValueError:      ## (skip tick if the ticker.MaxNLocator gave invalid target_value to brentq optimization)
        pass
## Finally, we set the positions and tick values
ax3.set_xticks(valid_xtick2loc)
from matplotlib.ticker import FixedFormatter
ax3.xaxis.set_major_formatter(FixedFormatter(["%g" % toptick for toptick in valid_xtick2val]))

## ==== Outputting ====
## Finish the plot + save 
ax.set_xlabel(u"Wavenumber $Ka/(2\pi)$");  # [m$^{-1}$]
ax3.set_xlabel(u"Wavenumber cm$^{-1}$");  # [m$^{-1}$]
ax.set_ylabel(u"Frequency $fa/(2\pi c)$") ## freq normalized 
ax2.set_ylabel(u"Frequency (THz)") ## freq normalized 
plt.legend(prop={'size':10}, loc='upper right')
plt.savefig("cdh_%s.pdf" % filesuffix, bbox_inches='tight')
plt.savefig("cdh_%s.png" % filesuffix, bbox_inches='tight')

    #plt.plot(freq, np.log10(np.abs(zf)+1e-10), color=color, label=u"$y'$", ls='-')      # (optional) plot amplitude
    #plt.plot(freq, np.unwrap(np.angle(zf)), color="#FF8800", label=u"$y'$", ls='--')   # (optional) plot phase

    ## FFT shift (to plot negative frequencies)

#if args.output[0:1] == '*':
    #outfilename = 'CDH_%s_%s' % (os.path.split(os.getcwd())[1], args.output[1:])
#else:
    #outfilename = args.output
#plt.savefig(outfilename, bbox_inches='tight')
