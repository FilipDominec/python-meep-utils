#!/usr/bin/env python
#coding:utf8    ## důležité: určuje kódování souboru; nutno zvolit cp1250 nebo utf8 dle editoru

import sys, gtk
import scipy,numpy
from numpy import pi
import matplotlib.pyplot as plt

## Make the Drude-Lorentz model
eps = 50. 
omega0 = 400e12           ## arbitrary low frequency that makes Lorentz model behave as Drude model
gamma = 100e12
omega_p = 2000e12        ## plasma frequency of aluminium
sigma = omega_p**2 / omega0**2

#{'omega': omega0, 'gamma': 25e12, 'sigma': omega_p**2 / omega0**2}, # (Lorentz) model for Al
def lorentz(x): #{{{
    return 1+sigma*omega0**2/(omega0**2 - 1j*x*gamma - x**2) 
#}}}
def stcMC(x, z, g, w):#{{{
    """ Stability check taken from MEEP comment
    
    Verified: if correct pole position given, returns zero at x=0+0j
    """
    A = (z + 1/z -2)
    B = g/2 * (z - 1/z)
    C = w**2
    return A*x**2 + B*x + C
    #x2 = (-B + numpy.sqrt(B**2 - 4*A*C)) / (2*A)
    #return abs(x2)
#}}}
def stcM(tsfreq, omega_0, gamma):#{{{
    """ Stability check from MEEP code 
    
    Behaves a bit weird, verify the stability on big computer
    """
    dt = 1/tsfreq/numpy.pi
    w = 2*numpy.pi*omega_0
    g = 2*numpy.pi*gamma;
    g2 = g*dt/2
    w2 = (w*dt)*(w*dt);
    b = (1 - w2/2) / (1 + g2)
    c = (1 - g2) / (1 + g2);
    return [b*b-c,    2*b*b - c + 2*abs(b)*numpy.sqrt(b*b - c) - 1] ## stable if one of them or both <0
#}}}
def z_magnitude(omega_0, gamma):#{{{
    """ Returns the the higher pole position in complex plane

    Lorentzian susceptibility is given as:
        sigma*omega0**2/(omega0**2 - 1j*x*gamma - x**2) 
    This function has two poles in the complex plane of x. Their position is
    determined by omega0 (resonance frequency) and gamma. 
    This function returns the position where the higher pole is found.
    """
    if omega_0 > gamma/2:
        return -1j*gamma/2 + numpy.sqrt(4*omega_0**2 - gamma**2) / 2
    else:
        return -1j*gamma/2 - 1j*numpy.sqrt(gamma**2 - 4*omega_0**2) / 2
#}}}
def print_clever_advice(omega_0, gamma, courant=.5):#{{{
    pass
#}}}

zm = z_magnitude(omega0, gamma)

print "ZM", zm, "magnitude", abs(zm), "max_res", 2.997e8/abs(zm)
re=numpy.arange(-15e14, 15e14, 1e13)#{{{
zmdata = []
zcdata = []
zdata = []
stc_data = []
imlist = []
for im in numpy.arange(-15e14, 15e14, 1e13):
    x = re+1j*im
    znew=lorentz(x)
    zmdata.append(abs(x)/abs(zm)-1)
    #zmdata.append(stcM(x, omega0, gamma)[0])
    #zmdata.append(stcM(x, omega0, gamma)[1])

    zcdata.append(stcMC(x, zm, 2*pi*omega0, 2*pi*gamma))

    zdata.append(znew)
    #stc_data.append(stc(x, omega0, gamma)[1])
    imlist.append(im)
zdata = numpy.array(zdata)
#stc_data = numpy.array(stc_data)
#}}}
plt.figure(figsize=(7,7))#{{{
## Plot eps polar
CS = plt.contourf(x, imlist, numpy.angle(zdata), numpy.arange(-numpy.pi+.001, numpy.pi, .01), cmap=plt.cm.jet, antialiased=True, lw=0)
CS2 = plt.contourf(x, imlist, numpy.arctan(abs(zdata)), numpy.arange(0., numpy.pi/2,  .01), cmap=plt.cm.gray, alpha=.5, antialiased=False, linewidths=0, lw=0)
#CS2 = plt.contour(x, imlist, abs(zdata),c='k')

## Plot eps cartesian
CS10 = plt.contour(x, imlist, numpy.real(zdata), numpy.arange(-10, 10), linewidths=0.5, colors='#0000aa', alpha=1)
plt.clabel(CS10, fontsize=7)
CS11 = plt.contour(x, imlist, numpy.imag(zdata), numpy.arange(-10, 10), linewidths=0.5, colors='#aa0000', alpha=1)
plt.clabel(CS11, fontsize=7)

## Plot the region of stability (should go through the pole position)
CS40 = plt.contour(x, imlist, zmdata, [0], colors='k', lw=2)
plt.clabel(CS40, fontsize=7)
plt.plot(zm.real, zm.imag, 'ks')

## Meep comment "... it is just a little algebra ..."
CS50 = plt.contour(x, imlist, zcdata, [0], colors='k', lw=3)
plt.clabel(CS50, fontsize=7)


#CS20 = plt.contour(x, imlist, abs(stc_data), 20, linewidths=0.5, colors='k', alpha=1)
#plt.clabel(CS20)
plt.xlabel(u"Re"); plt.ylabel(u"Im"); plt.legend(); plt.grid(); 
plt.savefig("cp.png")#}}}

## Plot the size of 

#re = numpy.arange(1e13,5e16,1e12)
plt.figure(figsize=(7,7))
znew=lorentz(re)
#plt.subplot(211)
#plt.subplot(212)
#plt.plot(re, stcMC(re, omega0, gamma), label="MEEP comment")
plt.plot(re, stcM(re, omega0, gamma)[0], label="MEEP code1")
plt.plot(re, stcM(re, omega0, gamma)[1], label="MEEP code2")
plt.plot(re, znew.real)
plt.plot(re, znew.imag, ls='--')
plt.plot(re, re-omega0, ls='--')
plt.ylim((-3,3)), plt.grid()
plt.legend()

plt.savefig("lp.png")
