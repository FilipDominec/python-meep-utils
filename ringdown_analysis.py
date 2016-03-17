#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
This is a manual computation of Fourier transform, to verify that the numpy.fft's or scipy's built-in Fast Fourier Transform
behaves as expected also in terms of accumulated power etc.

The curves obtained using the 
    1) frequency `f'
    2) angular frequency `omega'
are obviously different: when the _angular_ frequency `omega' is used, its axis dilates by 2pi, and 
accordingly, the values have to be divided by sqrt(2pi) (NB this is due to the Fourier-Plancherel theorem, 
which is tested below)

The non-angular frequency `f' approach is used by numpy.fft, and it can be shown that it gives data nearly identical to the manually computed FT.

Public domain, 2014 F. Dominec
"""

## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, hbar, pi

## == User settings ==
test_kramers_kronig = 0  # Compare data to their Hilbert transform (warning - short time record makes the KKR data ugly)
harmonic_inversion  = 1
analytic_lorentzian = 1  # Knowing the oscillator parametres, we can compare to the analytic solution
annotate            = 0  # Optionally add text labels

plot_absolute_value = 1
plot_ylog =           1

convention = 'f'            
#convention = 'omega'

FDMtrunc = (.0, 1.)         # Harminv (also known as FDM) may work better when it is supplied a shorter time record
                            # and it requires clipping the initial timespan when the broadband source is operating
                            # Note that if the time record is shorter than resonance lifetimes, FFT energy is less than that from Harminv

frequency_zoom = .05         # higher value = broader scope; use 1.0 to see all
frequency_full_axis = 0     # if enabled, does not clip the plotted range

# The following adds a Dirac-delta function at zero time, emulating e.g. the vacuum permittivty
# (note: Plancherel theorem is not valid for delta function in numerical computation)
add_delta_function =  0

## Prepare plot, optionally use LaTeX
#matplotlib.rc('text', usetex=True)
#matplotlib.rc('font', size=12)
#matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{txfonts}, \usepackage{lmodern},')
#matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})  ## select fonts
plt.figure(figsize=(10,5))


## == Time domain ==
plt.subplot(121)

if len(sys.argv) <= 1: 
    ## Generate time-domain data
    x, omega0, gamma, ampli = np.linspace(0., 10e-12, 4000), 2*pi*5e12, 2*pi*3e11, 1. ## note: everything should work for any frequency scale
    #x, omega0, gamma, ampli = np.linspace(0, 25, 3000), 2*pi*2, 2*pi*.3, 1.

    y = ampli * (np.sign(x)/2+.5) * np.sin(x*omega0) * np.exp(-x*gamma/2)           ## damped oscillator
    ## Note: Delta function is suitable for f-convention only (normalize by 2pi if omega-convention is used)
    if add_delta_function:
        y[int(len(x)*(-x[0]/(x[-1]-x[0])))] += 1 / (x[1]-x[0])  
    analytic_input = True
else:
    ## Load time-domain data
    try:
        data = np.loadtxt(sys.argv[1], unpack=True)
        if len(data) == 3:
            x, Eabs, Ephase = data
            y = Eabs * np.exp(1j*Ephase) # TODO harminv fails to find heavily damped oscillators;   to test out, add something like: * np.exp(-x/1e-12)
        else:
            x, y = data
        analytic_input = False
    except IndexError:
        print "Error: if a timedomain file is provided, it must have 2 or 3 columns: time, amplitude, [phase]"; quit()

maxplotf = frequency_zoom / (x[1]-x[0])

## Plot time-domain
plt.plot(x, y.real, c='#aa0088', label="Real part")
plt.plot(x, y.imag, c='#aa0088', label="Imaginary part", ls='--')
plt.grid(); plt.yscale('linear'); 
plt.legend(prop={'size':10}, loc='upper right'); plt.xlabel(u"time $t$"); plt.ylabel(u"response"); plt.title(u"a) Time domain")
Wt = np.trapz(y=np.abs(y)**2, x=x); print 'Plancherel theorem test: Energy in timedomain              :', Wt


## == Frequency domain ==
plt.subplot(122)

## An exact curve for the analytic solution of a damped oscillator
def lorentz(omega, omega0, gamma, ampli):
    return ampli / (omega0**2 - omega**2 + 1j*omega*gamma) 


def plot_complex(x, y, **kwargs):
    if plot_absolute_value:
        plt.plot(x, np.abs(y), **kwargs)
    else:
        kwargsr = kwargs.copy(); kwargsr['label']+=' (real)'; plt.plot(x, y.real, **kwargsr)
        kwargsi = kwargs.copy(); kwargsi['label']+=' (imag)'; plt.plot(x, y.imag, ls='--', **kwargsi)

## Spectrum 1: Scipy's  implementation of Fast Fourier transform
freq    = np.fft.fftfreq(len(x), d=(x[1]-x[0]))                 # calculate the frequency axis with proper spacing
yf2     = np.fft.fft(y, axis=0) * (x[1]-x[0])                   # calculate FFT values (maintaining the Plancherel theorem)
freq    = np.fft.fftshift(freq)                                 # reorders data to ensure the frequency axis is a growing function
yf2     = np.fft.fftshift(yf2) / np.exp(1j*2*pi*freq * x[0])    # dtto, and corrects the phase for the case when x[0] != 0
truncated = np.logical_and(freq>-maxplotf, freq<maxplotf)         # (optional) get the frequency range
(yf2, freq) = map(lambda x: x[truncated], (yf2, freq))    # (optional) truncate the data points
plot_complex(freq, yf2, c='#ff4400', label='ScipyFFT in $f$', lw=2, alpha=1)
Ws = np.trapz(y=np.abs(yf2)**2, x=freq); print 'Plancherel theorem test: Energy in freqdomain f (by Scipy) :', Ws 

## Spectrum 2: Own implementation of slow Fourier transform
f = np.linspace(-maxplotf, maxplotf, 1000)
yf = np.sum(y * np.exp(-1j*2*pi*np.outer(f,x)), axis=1) * (x[1]-x[0])
plot_complex(f, yf, c='#ffaa00', label='Manual FT in $f$', lw=.5, alpha=1)
Wm = np.trapz(y=np.abs(yf)**2, x=f); print 'Plancherel theorem test: Energy in freqdomain f (manual)   :', Wm

## Spectrum 3: Hilbert transform of the Fourier-transformed signal should be close to 1j times the signal
def naive_hilbert_transform(x, y, new_x): ## or, just a discrete convolution with the 1/t function
    old_x_grid, new_x_grid = np.meshgrid(x, new_x)
    sharpness = 5000         # with ideally dense sampled data, this should converge to infinity; reduce it to avoid ringing 
    return -1j * np.sum(y * np.arctan(1/(new_x_grid - old_x_grid)/sharpness)*sharpness, axis=1) / len(x) / (2*pi)
if test_kramers_kronig:
    ## Test the Kramers-Kronig relations - in f
    new_f = np.linspace(-maxplotf, maxplotf, 1000)
    conv = naive_hilbert_transform(f, yf, new_f)
    #conv = naive_hilbert_transform(freq, yf2, new_f) ## FIXME - using these data, the KK amplitude goes wrong
    plot_complex(new_f, conv, c='b', alpha=.6, lw=.5, label='KKR in $f$') 

## Spectrum 4: use the Filter-Diagonalisation Method, implemented in Harminv, to find discrete oscillators
if harmonic_inversion:
    import harminv_wrapper
    tscale = 1.0     ## harminv output may have to be tuned by changing this value
    x = x[int(len(x)*FDMtrunc[0]):int(len(x)*FDMtrunc[1])]*tscale
    y = y[int(len(y)*FDMtrunc[0]):int(len(y)*FDMtrunc[1])]
    hi = harminv_wrapper.harminv(x, y, amplitude_prescaling=None)
    hi['frequency'] *= tscale 
    hi['decay'] *= tscale

    oscillator_count = len(hi['frequency'])
    freq_fine = np.linspace(-maxplotf, maxplotf, 2000)
    sumosc = np.zeros_like(freq_fine)*1j
    for osc in range(oscillator_count):
        osc_y = lorentz(omega=freq_fine*2*pi,   
                omega0=hi['frequency'][osc]*2*pi, 
                gamma=hi['decay'][osc]*2*pi, 
                ampli=hi['amplitude'][osc] * np.abs(hi['frequency'][osc])*2*pi)
        sumosc += osc_y 
    plot_complex(freq_fine, sumosc, color='g', alpha=.6, lw=2, label=u"Harminv modes sum")      # (optional) plot amplitude
    Wh = np.trapz(y=np.abs(sumosc)**2, x=freq_fine); print 'Plancherel theorem test: Energy in Harminv f               :', Wh, '(i.e. %.5g of timedomain)' % (Wh/Wt)
    print 'All harminv oscillators (frequency, decay and amplitude):\n', np.vstack([hi['frequency'], hi['decay'], hi['amplitude']])

## Spectrum 5: If one Lorentzian is defined, we can plot its shape directly
if analytic_input and analytic_lorentzian:
    lor = lorentz(omega=f*2*pi, omega0=omega0, gamma=gamma, ampli=ampli*omega0)
    plot_complex(f, lor, c='b', alpha=.3, lw=5,  label='Exact Lorentzian in $f$') 
    Wa = np.trapz(y=np.abs(lor)**2, x=f); print 'Plancherel theorem test: Energy in analytic osc (f)        :', Wa 
    print 'Analytic    oscillators frequency, decay and amplitude:\n', np.vstack([omega0/2/pi, gamma/2/pi, ampli])

if annotate:
    try:
        with open('./annotate.txt') as f:
            annotations = {}
            for line in f:      text, freq = line.split('\t', 1); annotations[float(freq)] = text
            import meep_utils
            meep_utils.annotate_frequency_axis(annotations, label_position_y=np.sum(np.abs(yf[0:2000]))/len(yf)/3, arrow_length=3, log_y=True)
            #from meep_utils import annotate_frequency_axis #TODO
            #annotate_frequency_axis(annotations, label_position_y=np.sum(np.abs(yf[0:2000]))/len(yf)/3, arrow_length=3, log_y=True)
    except IOError: 
        print 'Error: the file ./annotate.txt could not be found'



## Finish the frequency-domain plot + save 
if not frequency_full_axis: plt.xlim(left=0, right=maxplotf)
plt.xscale('linear')
#plt.ylim((-16,16)); 
if plot_ylog:
    plt.yscale('log')
else:
    plt.yscale('linear')

if convention == 'f':
    plt.xlabel(u"frequency $f$"); 
elif convention == 'omega': 
    plt.xlabel(u"angular frequency $\\omega$"); 
plt.ylabel(u"local permittivity $\\varepsilon^{\\rm(Loc)}(\\omega)$"); 
plt.title(u"b) Frequency domain")
plt.grid()

plt.legend(prop={'size':10}, loc='upper right')
plt.savefig("oscillator_spectrum.png", bbox_inches='tight')

