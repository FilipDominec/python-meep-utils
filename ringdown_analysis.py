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

The frequency `f' approach is used by numpy.fft, and it can be shown that it gives data nearly identical to the manually computed FT.

Public domain, 2014 F. Dominec
"""

## Import common moduli
import matplotlib, sys, os, time
import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, hbar, pi

## == User settings ==

convention = 'f'            
#convention = 'omega'

# Compare data to their Hilbert transform (warning - short time record makes the KKR data ugly)
test_kramers_kronig = 0   

# The following adds a Dirac-delta function at zero time, emulating e.g. the vacuum permittivty
# (note: Plancherel theorem is not valid for delta function in numerical computation)
add_delta_function =  0

harmonic_inversion  = 1

# Knowing the oscillator parametres, we can compare to the analytic solution
analytic_lorentzian = 0


## == /User settings ==



## Use LaTeX
#matplotlib.rc('text', usetex=True)
#matplotlib.rc('font', size=12)
#matplotlib.rc('text.latex', preamble = '\usepackage{amsmath}, \usepackage{yfonts}, \usepackage{txfonts}, \usepackage{lmodern},')
#matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern Roman, Times']})  ## select fonts
plt.figure(figsize=(20,10))
plt.subplot(121)


## Generate time-domain data


x, omega0, gamma, ampli = np.linspace(-0e-3, 45e-3, 4000), 2*np.pi*1e3*1, 2*np.pi*.03*1e3, 1.
#x, omega0, gamma, ampli = np.linspace(-0, 2.5, 3000), 20*np.pi, 20*np.pi*.1, 1.
#x, omega0, gamma, ampli = np.linspace(-3, 25, 3000), 2*np.pi*3, 2*np.pi*.3, 1.
maxplotf = 200 / np.max(x)
y = ampli * (np.sign(x)/2+.5) * np.sin(x*omega0) * np.exp(-x*gamma/2)  ## damped oscillator
#y +=ampli * (np.sign(x)/2+.5) * np.sin(x*omega0*3)*2*pi * np.exp(-x*gamma/2)  ## damped oscillator
if add_delta_function:
    if convention == 'f':
        y[int(len(x)*(-x[0]/(x[-1]-x[0])))] +=         1 / (x[1]-x[0])  ## delta function suitable for f-convention 
    elif convention == 'omega':
        print "Warning: Delta function unclear how to be implemented in omega convention"
        y[int(len(x)*(-x[0]/(x[-1]-x[0])))] +=         1 / (x[1]-x[0])  ## delta function suitable for omega-convention 


## Plot time-domain
plt.plot(x,y, c='#aa0088', label="Real part")
plt.grid()
plt.yscale('linear')

## Prepare the analytic solution
def lorentz(omega, omega0, gamma, ampli):
    return ampli / (omega0**2 - omega**2 + 1j*omega*gamma) 
    #return ampli * omega0**2 / (omega0**2 - omega**2 + 1j*omega*gamma) 

## Plot time-domain 
plt.xlabel(u"time $t$") 
plt.ylabel(u"medium response function $\\chi_e^{\\rm(Loc)}(t) + \\delta(t)$") 
plt.title(u"\\textbf{a)} Time domain")
print 'Plancherel theorem test: Energy in timedomain              :', np.trapz(y=np.abs(y)**2, x=x)

def naive_hilbert_transform(x, y, new_x):
    old_x_grid, new_x_grid = np.meshgrid(x, new_x)
    return -1j * np.sum(y * np.arctan(1/(new_x_grid - old_x_grid)/200)*200, axis=1) / len(x) * 2*pi

plt.subplot(122)
if convention == 'f':
    ## Scipy's  implementation of Fast Fourier transform
    freq    = np.fft.fftfreq(len(x), d=(x[1]-x[0]))                 # calculate the frequency axis with proper spacing
    yf2     = np.fft.fft(y, axis=0) * (x[1]-x[0])                   # calculate FFT values (maintaining the Plancherel theorem)
    freq    = np.fft.fftshift(freq)                                 # reorders data to ensure the frequency axis is a growing function
    yf2     = np.fft.fftshift(yf2) / np.exp(1j*2*pi*freq * x[0])    # dtto, and corrects the phase for the case when x[0] != 0
    truncated = np.logical_and(freq>-maxplotf, freq<maxplotf)         # (optional) get the frequency range
    (yf2, freq) = map(lambda x: x[truncated], (yf2, freq))    # (optional) truncate the data points
    plt.plot(freq, yf2.real, c='m', label='ScipyFT in $f$', lw=2)
    plt.plot(freq, yf2.imag, c='m', ls='--', lw=2)
    Ws = np.trapz(y=np.abs(yf2)**2, x=freq); print 'Plancherel theorem test: Energy in freqdomain f (by Scipy) :', Ws 

    ## Own implementation of slow Fourier transform - in f
    f = np.linspace(-maxplotf, maxplotf, 1000)
    yf = np.sum(y * np.exp(-1j*2*pi*np.outer(f,x)), axis=1) * (x[1]-x[0])
    plt.plot(f, yf.real, c='g', label='ManFT in $f$')
    plt.plot(f, yf.imag, c='g', ls='--')
    Wm = np.trapz(y=np.abs(yf)**2, x=f); print 'Plancherel theorem test: Energy in freqdomain f (manual)   :', Wm

    if test_kramers_kronig:
        ## Test the Kramers-Kronig relations - in f
        new_f = np.linspace(0, 5, 100)
        conv = naive_hilbert_transform(f, yf, new_f)
        plt.plot(new_f, conv.real, ls='-', c='k', alpha=1, lw=.5, label='KKR in $f$') 
        plt.plot(new_f, conv.imag, ls='--', c='k', alpha=1, lw=.5) 

    if analytic_lorentzian:
        lor = lorentz(omega=f*2*pi, omega0=omega0, gamma=gamma, ampli=ampli*omega0)
        plt.plot(f, lor.real, ls='-',  c='r', alpha=.5, lw=1.5,  label='Osc in $f$') 
        plt.plot(f, lor.imag, ls='--', c='r', alpha=.5, lw=1.5) 
        Wa = np.trapz(y=np.abs(lor)**2, x=f); print 'Plancherel theorem test: Energy in analytic osc (f)        :', Wa 

    if harmonic_inversion:
        import harminv_wrapper
        amplitude_prescaling = 1e3
        hi = harminv_wrapper.harminv(x, y, amplitude_prescaling=amplitude_prescaling)
        #hi = harminv_wrapper.harminv(x[0:int(len(x)*.6)], y[0:int(len(x)*.6)], amplitude_prescaling=3.)
        if type(hi['frequency']) == np.float64:
            hi['frequency'], hi['decay'], hi['amplitude'] = np.array([hi['frequency']]), np.array([hi['decay']]), np.array([hi['amplitude']])
        print hi['frequency'], hi['decay'], hi['amplitude']
        oscillator_count = len(hi['frequency'])
        freq_fine = np.linspace(-maxplotf, maxplotf, 2000)
        sumosc = np.zeros_like(freq_fine)*1j
        for osc in range(oscillator_count):
            #osc_y = lorentz(omega=freq_fine*2*pi,   omega0=hi['frequency'][osc]*2*pi, gamma=hi['decay'][osc]*4, ampli=hi['amplitude'][osc]*pi**2)
            osc_y = lorentz(omega=freq_fine*2*pi,   
                    omega0=hi['frequency'][osc]*2*pi, 
                    gamma=hi['decay'][osc]*4, 
                    ampli=hi['amplitude'][osc] * np.abs(hi['frequency'][osc]) * np.abs(hi['decay'][osc]) / 39.**.5) 
            sumosc += osc_y 
        plt.plot(freq_fine, sumosc.real, color="#0088FF", label=u"$\\Sigma$ osc", ls='-')      # (optional) plot amplitude
        plt.plot(freq_fine, sumosc.imag, color="#00aaFF", label=u"", ls='--')      # (optional) plot amplitude
        Wh = np.trapz(y=np.abs(sumosc)**2, x=freq_fine); print 'Plancherel theorem test: Energy in Harminv f   :', Wh 
        print "Wh/Ws", Wh/Ws

elif convention == 'omega':
    # Own implementation of slow Fourier transform - in omega XXX
    omega = np.linspace(-maxplotf*2*pi, maxplotf*2*pi, 3000)  # (note: if only positive frequencies are used, the energy will be half of that in time-domain)
    yomega = np.sum(y * np.exp(-1j*        np.outer(omega,x)), axis=1) * (x[1]-x[0])  / np.sqrt(2*pi)
    plt.plot(omega, np.real(yomega), c='#440088', label='Real part') # , label='FT in $\\omega$-convention'
    plt.plot(omega, np.imag(yomega), c='#0088aa', lw=.8, ls='--', label='Imaginary part')
    print 'Plancherel theorem test: Energy in freqdomain omega :', np.trapz(y=np.abs(yomega)**2, x=omega)

    if test_kramers_kronig:
        ## Test the Kramers-Kronig relations - in omega
        new_omega = np.linspace(5, 8, 500)
        conv = naive_hilbert_transform(omega, yomega, new_omega)
        plt.plot(new_omega, conv.real, ls='-', c='k', alpha=1, lw=.5, label='KKR in $f$', ms=3, marker='o') 
        plt.plot(new_omega, conv.imag, ls='--', c='k', alpha=1, lw=.5) 

    if analytic_lorentzian:
        lor = lorentz(omega=omega, omega0=omega0, gamma=gamma, ampli=ampli) / (2*pi)**.5
        plt.plot(omega, lor.real, ls='-',  c='r', alpha=1, lw=.5,  label='Osc in $f$') 
        plt.plot(omega, lor.imag, ls='--', c='r', alpha=1, lw=.5) 

## Finish the frequency-domain plot + save 
#plt.xlim(left=0)
plt.xscale('linear')
#plt.ylim((-16,16)); 
plt.yscale('linear')
if convention == 'f':
    plt.xlabel(u"frequency $f$"); 
elif convention == 'omega': 
    plt.xlabel(u"angular frequency $\\omega$"); 
plt.ylabel(u"local permittivity $\\varepsilon^{\\rm(Loc)}(\\omega)$"); 
plt.title(u"\\textbf{b)} Frequency domain")
plt.grid()

plt.legend(prop={'size':10}, loc='upper right')
plt.savefig("oscillator_spectrum.png", bbox_inches='tight')

