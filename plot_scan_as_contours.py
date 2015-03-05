#!/usr/bin/env python
#-*- coding: utf-8 -*-

## Import common moduli
import sys, os.path
import numpy as np
from scipy.constants import c, hbar, pi
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm ## XXX TODO
import meep_utils 


## ==== User settings ====

## Select which structure properties will be plot and which not (each in separate file)
quantities = []
quantities += ['ampli']

parameter_name = 'height'
quantity = 'amplitude'
xlabel = 'Frequency (GHz)'
ylabel = 'Cylinder shape $(2a/h)^2$'
recalculate_to_angle = False
logarithmic = 0
frequnit = 1e9
yunit =  1
minf, maxf  = 0, 7e9           # span of the horizontal axis to be plot
xlim        = (minf/frequnit, maxf/frequnit)
ylim = (0, .02)
interp_anisotropy = 1       # value lower than 1. interpolates rather vertically; optimize if plot desintegrates

def y_function(y):
    return (2*parameters['radius']/y)**2 
usetex = 0
## ==== / User settings ====

if usetex:
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=12)
    matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{upgreek}, \usepackage{palatino}')
    matplotlib.rc('font',**{'family':'serif','serif':['palatino, times']})  ## select fonts   , 
    #matplotlib.rc('font',**{'family':'sans-serif', 'sans-serif':['Computer Modern Sans serif']})

##Start figure + subplot
plt.figure(figsize=(15,11))

## Load data from multiple files
if len(sys.argv) > 1:
    filenames = sys.argv[1:]
else: 
    filenames = [x for x in os.listdir(os.getcwd())   if '.dat' in x]
print "Got %d files to plot %s..." % (len(filenames), quantity)
x, y, z = [np.array([]) for _ in range(3)] ## three empty arrays

#filenames.sort(key=lambda name: float(name.split('radius=')[1].split('_')[0]))     # sort (optional)

for datafile_name in filenames: 
    ## Load header
    parameters  = meep_utils.loadtxt_params(datafile_name)
    columns     = meep_utils.loadtxt_columns(datafile_name)

    ## Getting 1D data
    try:
        (freq, ampli, phase) = np.loadtxt(datafile_name, usecols=range(3), unpack=True)
        #(freq, s11_ampli, s11p, s12_ampli, s12p, Nre, Nim, Zre, Zim, eps_r, eps_i, mu_r, mu_i) = \
                #np.loadtxt(datafile_name, usecols=range(13), unpack=True)
    except IndexError:
        print 'Error: Insufficient number of columns in data file'; quit()

    ## Convert the data to a position in the plot, optionally making some transformations
    input_y = parameters['height']
    input_z = ampli
    znew = np.log10(input_z)
    ynew = y_function(input_y) * np.ones_like(freq) 

    ## Store the points
    x = np.append(x, freq/frequnit)
    y = np.append(y, ynew/yunit)
    z = np.append(z, znew)

if not ylim: ylim=(min(y), max(y))
plt.ylabel(ylabel); 
plt.xlabel(xlabel);

# Grid the data.
from matplotlib.mlab import griddata
xi = np.linspace(min(x), max(x), 600)
yi = np.linspace(min(y), max(y), 500)
zi = griddata(x, y*interp_anisotropy, z, xi, yi*interp_anisotropy, interp='linear')

# Standard contour plot
levels = np.linspace(np.min(zi), np.max(zi), 100)
cmap = cm.gist_earth;
contours = plt.contourf(xi,yi,zi, cmap=cmap, levels=levels, extend='both')  
for contour in contours.collections: contour.set_antialiased(False) ## fix aliasing for old Matplotlib
plt.colorbar().set_ticks(list(range(0, int(np.max(levels)+1))))

## Plot analytic mode frequencies ## XXX for cylindrical resonator only
#from scipy.special import jnyn_zeros
#freq_correction = (1. - 350e6/3.2e9)
#radius, height = parameters['radius'], parameters['height']
#plot_height = np.linspace(120e-3, 1000e-3, 100)
#for p in [0,1,2,3,4]:
    #for n in range(4):
        #S = " "*p
        #for m,B in enumerate(jnyn_zeros(n, 5)[0]):
            #analytic_freq = freq_correction * c/(2*np.pi) * np.sqrt((B/(parameters['radius']))**2 + (p*np.pi/plot_height)**2)
            ## if analytic_freq/frequnit < 10:
            #plt.plot(analytic_freq/frequnit, y_function(plot_height)/yunit, c='b', lw=.8)
        #for m,B in enumerate(jnyn_zeros(n, 5)[1]):
            #if p>0: ## TExx0 modes can not exist [Pozar]
                #analytic_freq = freq_correction * c/(2*np.pi) * np.sqrt((B/(parameters['radius']))**2 + (p*np.pi/plot_height)**2)
            ## if analytic_freq/frequnit < 10:
            #plt.plot(analytic_freq/frequnit, y_function(plot_height)/yunit, c='k', lw=.8)

plt.xlim(xlim)
plt.grid()
plt.savefig('%s_%s.png' % (quantity, os.path.split(os.path.dirname(os.getcwd()))[1]), bbox_inches='tight')
