#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" Plots reflection and transmission of a metamaterial structure

Tries to calculate its effective parameters [Smith2002], avoiding branch jumps
Enables to save to several types of output (cartesian graphs, polar graphs, nice PDF graphs...)
Exports the effective parameters to another data file for further processing

About this script:
 * Written in 2012-13 by Filip Dominec (dominecf at the server of fzu.cz).
 * Being distributed under the GPL license, this script is free as speech after five beers. 
 * You are encouraged to use and modify it as you need. Feel free to write me if needed.
 * Hereby I thank to the MEEP/python_meep authors and people of meep mailing list who helped me a lot.

TODOs:
 * Guess the correct branch for N (using Kramers-Kronig relations?)
 * Fix the passivity criterion for Im N > 0, Re Z > 0

"""
import numpy as np
import sys, os, re, matplotlib 
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, fmin
from scipy.constants import pi, c

## == User settings for postprocessing and plotting == 
frequnit, frequnitname = 1e12, "THz"

N_init_branch   =   -1
N_init_sign     =   1
autocorrect_signs = True
Z_init_sign     =   -1
check_hilbert   =   0       ## Verifies if Kramers-Kronig relations hold for N  ###XXX
legend_enable   =   1      
brillouin_boundaries = 1    ## Plots thin lines where the N would exceed the allowed 
                            ## range for 0-th Bloch mode
autobranch      = 0

savedat     = 1     ## created directory 'effparam' and saves all params to an ascii file with header
plot_publi  = 0     ## prepares nice small graphs for publication
plot_polar  = 0     ## plots results to polar graphs for diagnostics
plot_bands  = 0     ## plots index of refraction as dispersion curves (k-omega)
plot_expe   = 1     ## if 'r.dat', 't.dat', 'N.dat', 'Z.dat', 'eps.dat' or 'mu.dat' available, overlay them
find_plasma_frequency = 0 ## find frequencies where epsilon crosses zero

plot_freq_min = 0
plot_freq_max = None  ## if None, decide from the input file header
#plot_freq_max = 2.5e12
padding = None
autobranch_sampler_position = 0.06

np.seterr(all='ignore')      ## do not print warnings for negative-number logarithms etc.
## == </user settings> == 

## == Auxiliary functions ==
def get_simulation_name(argindex=1): #{{{
    """Get the name of the last simulation run.

    Priority: 1) parameter, 2) last_simulation_name.dat, 3) working directory"""
    cwd = os.getcwd()
    if len(sys.argv)>argindex and sys.argv[argindex] != "-"  and __name__ == "__main__": 
        print "Parameter passed:", sys.argv[argindex]
        last_simulation_name = sys.argv[argindex]
    elif os.path.exists(os.path.join(cwd, 'last_simulation_name.dat')):
        print "Loading from", os.path.join(cwd, 'last_simulation_name.dat')
        last_simulation_name = os.path.join(cwd, open(os.path.join(cwd, 'last_simulation_name.dat'),'r').read().strip())
    else:
        print "Error: No input file provided and 'last_simulation_name.dat' not found!"
        last_simulation_name = cwd
    if (last_simulation_name[-4:] == ".dat"): last_simulation_name = last_simulation_name[:-4] # strip the .dat extension
    return  last_simulation_name
#}}}
def get_cmdline_parameters():#{{{ (unused)
    # (optional) Manual N branch override
    if len(sys.argv)>2 and sys.argv[2] != "-"  and __name__ == "__main__": 
        print "Setting branch:", sys.argv[2]
        branch_offset = np.ones(len(freq))*int(sys.argv[2])
        last_simulation_name += "_BRANCH=%s" % sys.argv[2]
    if len(sys.argv)>3 and sys.argv[3] != "-"  and __name__ == "__main__": 
        print "Setting branch sign:", sys.argv[3]
        Nsign = np.ones(len(freq))*int(sys.argv[3])
        last_simulation_name += "_SIGN=%s" % sys.argv[3]
    return branch_offset, Nsign#}}}
def load_rt(filename, layer_thickness=None, plot_freq_min=None, plot_freq_max=None, truncate=True, padding=None): #{{{
    """ Loads the reflection and transmission spectra and simulation settings 

    Returns:
    * frequency axis
    * reflection s11 and transmission s12 as complex np arrays

    Compatible with the PKGraph text data file with polar data: 
    * parameters in header like: #param name,value
    * column identification like: #column Ydata
    * data columns in ascii separated by space
    Expects polar data with columns: frequency, s11 ampli, s11 phase, s12 ampli, s12 phase
    """
    ## Extract relevant parameters
    with open(filename+'.dat') as datafile:
        for line in datafile:
            if line[0:1] in "0123456789": break         # end of file header
            value = line.replace(",", " ").split()[-1]  # the value of the parameter will be separated by space or comma
            if ("cell_size" in line) and (layer_thickness == None): cell_size = float(value)
            if ("cells" in line) and (layer_thickness == None): cells = float(value)
            if ("plot_freq_min" in line) and (plot_freq_min == None): plot_freq_min = float(value)
            if ("plot_freq_max" in line) and (plot_freq_max == None): plot_freq_max = float(value)
            if ("param padding" in line) and (padding == None): padding = float(value)

    ## Load data columns
    (freq, s11amp, s11phase, s12amp, s12phase) = \
            map(lambda a: np.array(a, ndmin=1), np.loadtxt(filename+".dat", unpack=True)) 
    
    ## If not specified, guess the plot frequency range
    if plot_freq_min == None: plot_freq_min = np.min(freq)
    if plot_freq_max == None: plot_freq_max = np.max(freq)

    ## Limit the frequency range to what will be plotted (recommended) TODO wrong approach
    if truncate and len(freq)>1:
        (d0,d1) = np.interp((plot_freq_min, plot_freq_max), freq, range(len(freq)))
        (freq, s11amp, s11phase, s12amp, s12phase) = \
                map(lambda a: a[int(d0):int(d1)], (freq, s11amp, s11phase, s12amp, s12phase))
    return freq, s11amp, s11phase, s12amp, s12phase, cell_size, plot_freq_min, plot_freq_max, padding, cells
#}}}
def shiftmp(freq, s11, shiftplanes):#{{{
    """ Adjusts the reflection phase like if the monitor planes were not centered.

    For symmetric metamaterial cell, this function is not needed. The symmetry requires that
    the monitor planes in front of and behind the mm cell are centered.

    However, for an asymmetric metamaterial, the correct position has to be found. Otherwise
    the Fresnel inversion gives negative imaginary part of N and/or negative real part of Z, which
    is quite spooky for passive medium. 
    
    Even such metamaterials, however, may be properly homogenized if we define the 
    position of monitor planes as a function of frequency. We can assume that:
    1) This optimum shift shall hold for all simulations with one or more unit cells.
    2) When the wave direction is reversed (to get s21, s22 parameters), the shift should be negated.
    These rules should enable us to homogenize any asymmetric non-chiral metamaterial. 

    Note that this shifting is still an experimental technique and has to be tested out thoroughly. 
    """
    return np.array(s11) * np.exp(1j*np.array(shiftplanes)/(c/freq) * 2*pi * 2)
#}}}
def find_maxima(x, y, minimum_value=.1):#{{{
    """ 
    Returns the x points where 
    1) y has a local maximum (i. e. dx/dy goes negative) AND 
    2) where y is above minimum_value 
    """
    d = y[1:-1] - y[0:-2]   ## naïve first derivative
    maxima = x[1:][np.sign(d[0:-2])-np.sign(d[1:-1]) + np.sign(y[2:-2]-minimum_value)==3]
    return maxima 
#}}}
def reasonable_ticks(a, density=.6): #{{{
    """ Define the grid and ticks a bit denser than by default """
    decimal=10**np.trunc(np.log10(a/density)); y=a/density/decimal/10
    return (decimal, 2*decimal, 5*decimal)[np.int(3*y)]
#}}}

## == Homogenisation functions (processing whole np.array at once) ==
def polar2complex(amp, phase): return amp*np.exp(1j*phase)  #{{{
#}}}
def unwrap_ofs(p, ofs):#{{{
    """ Similar to np.unwrap, but take into account the initial offset. 
    Increment this offset if needed, and return it as the second return value.
    """
    return np.unwrap(p)+ofs, (np.unwrap(p)-p)[-1]+ofs
#}}}
def rt2n(frequency, s11, s12, d, init_branch=0, init_sign=1, uo=[0,0,0,0]): #{{{
    """ Invert Fresnel equations to obtain complex refractive index N, with autodetection of arccosine branch#{{{

    Accepts:
    * s11 - np.array of reflection,
    * s12 - np.array of transmission,
    * d   - layer thickness

    Returns: a tuple of three np.arrays
      * the retrieved effective index of refraction, 
      * the arccos() branch used for its calculation, 
      * the debug information 
    Technical details are commented in the code.
    
    This algorithm implements the method for effective refractive index retrieval from the 
    s11 and s12 scattering coefficients [Smith2002]. 

    Such calculation is not unambiguous due to multiple branches of complex arccosine. If the branches
    of the solution are not switched at proper frequencies, the index of refraction often becomes discontinuous
    and generally wrong.

    This function fixes this problem by the analysis of the arccos() behaviour. It requires that the r(f) and t(f)
    are supplied as whole spectra. It is then possible to trace the argument passed to arccos() and properly choose
    the correct branch for whole spectral region.


    Limitations of this procedure:
        * If structure is highly reflective at lowest frequencies (metallic wires etc.), the N branch cannot be determined
          reliably. To fix this, increase 'plot_freq_min' (the start of computed frequency range), or provide init_branch.

          Initial branch choosing is not implemented. Its value may be optionally provided in the argument init_branch and 
          init_sign. The user should choose theme so that the initial values for 
            i) the curves are continuous 
            ii) and: Im(N) > 0 (for a nonamplifying medium)

        * The np.unwrap() function requires that the frequency is sampled fine enough. If the branch is wrongly detected
          at sharp resonances, there are good news: You probably do not have to run the simulation longer; often is 
          sufficient to pad the time-domain data with zeroes.

        * For some simulations, there is a weird _continuous_ branch transition at higher frequencies for thicker 
          metamaterial samples. The returned index of refraction breaks Kramers-Kronig relations.
          However, the Hilbert transform of the imaginary part of N" gives proper data. Single-cell simulation also gives 
          proper data...

          Putting the layers far apart alleviates this for 2 cells: can it be related to higher-order Bloch modes? 

    """#}}}

    ## Argument passed to arccos():
    arg = (1+0j-s11**2+s12**2)/2/(s12)

    ## Count passing through complex arccos() branch cuts in the complex plane:
    lu, uo[0] = unwrap_ofs(np.angle(arg + 1. + 0e-3j) + pi, uo[0])
    ru, uo[1] = unwrap_ofs(np.angle(arg - 1. + 0e-3j), uo[1])

    lbc = np.floor(lu/2/pi)            
    rbc = np.floor(ru/2/pi)
    anl = (-1)**(lbc) ## left cut:  (-inf .. -1]
    anr = (-1)**(rbc) ## right cut: [1 .. +inf)

    ## Retrieve the sign and branch of the arccos()
    sign = anr*anl*init_sign

    lbr, uo[2] = unwrap_ofs(np.angle(-anr + 1j*anl) + pi, uo[2])
    rbr, uo[3] = unwrap_ofs(np.angle(+anr - 1j*anl) + pi, uo[3])
    branch = np.floor(lbr/2/pi) + np.floor(rbr/2/pi) + 1 + init_branch
    #branch = np.floor(np.unwrap(np.angle(rbc + 1j*lbc))/2/pi) + \
            #np.floor(np.unwrap(np.angle(-rbc - 1j*lbc))/2/pi) + 1 + init_branch

    ## Standard Fresnel inversion:
    k = 2*pi * frequency/c          # the wave vector
    N = np.conj((np.arccos(arg)*sign + 2*pi*branch) / (k*d)) 

    #if abs(frequency[-1]-387.3e9)<1e9: ## debug
        #print "f branch uo", frequency, branch, uo
    return N, uo, (branch, sign, arg, anr, anl)
    """ For diagnostics, you may also wish to plot these values:#{{{

    #argLog = np.e**(1j*np.angle(arg))*np.log(1+abs(arg)) ## shrinked graph to see the topology

    plt.plot(freq, arg.real, color="#aaaaaa", label=u"$arg$'", lw=1) 
    plt.plot(freq, arg.imag, color="#aaaaaa", label=u"$arg$'", lw=1, ls='--') 
    #plt.plot(freq, argLog.real, color="#000000", label=u"$arg$'", lw=1) 
    #plt.plot(freq, argLog.imag, color="#000000", label=u"$arg$'", lw=1, ls="--") 
    #plt.plot(freq, np.ones_like(freq)*np.log(2),  color="#bbbbbb", label=u"$arg$'", lw=1) 
    #plt.plot(freq, -np.ones_like(freq)*np.log(2), color="#bbbbbb", label=u"$arg$'", lw=1) 
    #plt.plot(freq, anr, color="#aaaaff", label=u"$anr$'", lw=1) 
    #plt.plot(freq, anl, color="#aaffaa", label=u"$anr$'", lw=1) 
    #plt.plot(freq, anr_trunc, color="#0000ff", label=u"$anrR$'", lw=1) 
    #plt.plot(freq, anl_trunc*.9, color="#00dd00", label=u"$anrR$'", lw=1) 
    #plt.plot(freq, branch*.8, color="#dd0000", label=u"$turn$'", lw=2) 
    #plt.plot(freq, sign*.7, color="#ffbb00", label=u"$sign$'", lw=2) 
    """#}}}
#}}}
def rt2z(s11, s12, init_sign=1, uo=0):#{{{
    """ Invert Fresnel equations to obtain complex impedance Z

    This function complements the refractive index obtained by rt2n() with the effective impedance.

    The computation is much easier, because the only unambiguous function is the complex square root.
    It allows two solutions differing by their sign. To prevent discontinuities, we calculate the 
    square root in polar notation. 

    Initial sign may be supplied by the user.

    Returns complex impedance as np.array
    """

    #def get_phase(complex_data):
        #""" Unwraps and shifts the phase from Fourier transformation """
        #if len(complex_data) <= 1: return np.angle(complex_data)
        #phase, uo = unwrap,ofs(np.angle(complex_data), uo)
        #center_phase = phase[min(5, len(phase)-1)] ## 5 is chosen to avoid zero freq.
        #return phase-(round(center_phase/2/pi)*2*pi)

    ## Calculate square root arguments 
    Zarg1=((1+s11)**2 - s12**2)
    Zarg2=((1-s11)**2 - s12**2)

    ## Get impedance from polar notation of (Zarg1/Zarg2)
    Zamp   =  abs(Zarg1 / Zarg2)**.5           ## amplitude of square root
    if hasattr(Zarg1, '__len__') and len(Zarg1)>1:
        Zphase, uo = unwrap_ofs(np.angle(Zarg1/Zarg2), uo)      ## phase of square root (without discontinuities) TODO
    else:
        Zphase = np.angle(Zarg1/Zarg2)
        uo = 0
    Z = np.conj(np.exp(1j*Zphase/2) * Zamp) * init_sign
    return Z, uo

    """
    ### Possible experimental improvements:
    EnforceZrePos =     True 
    FlipZByPhaseMagic = True
    Zrealflipper = 1        ## unphysical if not 1
    Zconjugator =  1         

    ## Exception to the Re(Z)>0 rule:
    Z_turnaround = (-1)**np.round(Zphase/pi)
    if FlipZByPhaseMagic: 
        Z = Z * Z_turnaround
    ## For testing only
    Z = (Z.real * Zrealflipper + 1j*Z.imag * Zconjugator)  
    if EnforceZrePos:
        Z *= np.sign(Z.real)
    """
#}}}
def nz2epsmu(N, Z):#{{{
    """ Accepts index of refraction and impedance, returns effective permittivity and permeability"""
    return N/Z, N*Z
#}}}
def epsmu2nz(eps, mu):#{{{
    """ Accepts permittivity and permeability, returns effective index of refraction and impedance"""
    N = np.sqrt(eps*mu)
    N *= np.sign(N.imag)
    Z = np.sqrt(mu / eps)
    return N, Z
#}}}
def nz2rt(freq, N, Z, d):#{{{
    """ Returns the complex reflection and transmission parameters for a metamaterial slab.

    Useful for reverse calculation of eps and mu (to check results)
    
    Accepts: 
    * frequency array,
    * effective refractive index N, 
    * effective impedance Z, 
    * vacuum wave vector k and 
    * thickness d of the layer.
    """
    ## Direct derivation from infinite sum of internal reflections
    k = 2*pi * freq/c          # the wave vector
    t1 = 2 / (1+Z)              # transmission of front interface
    t2 = 2*Z / (Z+1)            # transmission of back interface
    t1prime = Z*t1       
    r1=(Z-1)/(Z+1)              # reflection of front interface
    r2=(1-Z)/(1+Z)              # reflection of back interface
    s12 = t1*t2*np.exp(1j*k*N*d) / (1 + r1*r2*np.exp(2j*k*N*d))
    s11 = r1 + t1prime*t1*r2*np.exp(2j*k*N*d)/(1+r1*r2*np.exp(2j*k*N*d))
    return s11, s12
    """
    Note: these results may be also re-expressed using goniometric functions.
    Equations from Smith2002 or Cai-Shalaev, mathematically equivalent to those above 
    (only Smith's s11 has negative sign convention).
    s12new = 1/(np.cos(N*k*d) - .5j*(Z+1/Z)*np.sin(N*k*d)) 
    s11new = -s12new * .5j*(Z-1/Z)*np.sin(N*k*d)      

    TODO: implement also for other surrounding media than vacuum.
    """
#}}}

## == Auxiliary functions for monitor-plane fitting == 
def error_func(N1,Z1,N2,Z2,lastdif=0,p0=[0]):#{{{
    """ Used for optimization: tries to match N1,N2 and Z1,Z2, avoiding forbidden values """
    return abs(N1-N2) + abs(Z1-Z2) + \
                lastdif + (abs(p0[0])*1e4)**2 +\
                (abs(np.imag(N1))-np.imag(N1))*100 + (abs(np.imag(N2))-np.imag(N2))*100 + \
                (abs(np.real(Z1))-np.real(Z1))*100 + (abs(np.real(Z2))-np.real(Z2))*100
                #}}}
def eval_point(p0):#{{{
    freq_p = freq[i-1:i+1]

    s11p1   = shiftmp(freq[i-1:i+1], s11[i-1:i+1], p0[0])
    s12p1   = s11[i-1:i+1] 
    new_N1, Nuo1x = rt2n(freq_p, s11p1, s12p1, d, init_branch=0, uo=Nuo1)[0:2]
    new_Z1, Zuo1x = rt2z(s11p1, s12[i-1:i+1], uo=Zuo1)

    s11p2   = shiftmp(freq[i-1:i+1], s11_2[i-1:i+1], p0[0])
    s12p2   = s11_2[i-1:i+1] 
    new_N2, Nuo2x = rt2n(freq_p, s11p2, s12p2, d2, init_branch=0, uo=Nuo2)[0:2]
    new_Z2, Zuo2x = rt2z(s11p2, s12[i-1:i+1], uo=Zuo2)

    lastdif = abs(p0s[-1]-p0[0])*1e5 if (p0s[-1] != np.NaN) else 0
    return error_func(new_N1[1], new_Z1[1], new_N2[1], new_Z2[1], lastdif=lastdif)
#}}}


## --- Calculation -------------------------------------------- 
## Get reflection and transmission data
last_simulation_name = get_simulation_name()
freq, s11amp, s11phase, s12amp, s12phase, cell_size, plot_freq_min, plot_freq_max, padding, cells = \
        load_rt(last_simulation_name, plot_freq_min=plot_freq_min, plot_freq_max=plot_freq_max, truncate=False, padding=padding)

d = cell_size * cells

## Convert to complex numbers and compensate for the additional padding of the monitor planes
s11 = shiftmp(freq, polar2complex(s11amp, s11phase), padding*np.ones_like(freq))
s12 = shiftmp(freq, polar2complex(s12amp, s12phase), padding*np.ones_like(freq))

## Build the debug plots
arg = (1+0j-s11**2+s12**2)/2/(s12)
argLog = np.e**(1j*np.angle(arg))*np.log(1+abs(arg)) ## radially shrinked graph to see the topology

## Calculate N, Z and try to correct the signs (TODO use K-K branch selection!)
if len(freq)>2:
    N, N_uo, N_debug    = rt2n(freq, s11, s12, d, init_branch=N_init_branch, init_sign=N_init_sign)
    #print "N before correctio1", N[0:10]
    Z, Z_uo             = rt2z(s11, s12, init_sign=Z_init_sign)
    if autocorrect_signs: 

        ## Fix N sign so that N.imag > 0 
        nimag = np.sum(np.clip(N.imag,-10., 10.))
        if nimag<0: N *= -1

        ## Take a sample
        ii = int(float(len(N)) * autobranch_sampler_position)
        sampleR = np.real(N[ii])
        sampleI = np.imag(N[ii])
        print 'sampleR, sampleI =', sampleR, sampleI

        #if ii == -0: N *= -1  ## alternative way of fixing sign?

        ## Fix N branch so that N.real does not diverge  at low frequencies
        branch_selector = 2*sampleR*freq[ii]/c*d
        det_branch = np.round(branch_selector)
        if (branch_selector-det_branch) < 0:
            N *= -1
            det_branch *= -1
            branch_selector *= -1
        print 'branch_selector, det_branch, diff', branch_selector, det_branch, (branch_selector-det_branch)
        #print "N before correction", N[0:20]
        N -= det_branch / (freq/c*d)/2
        #print "N after  correction", N[0:20]
        ## Fixing Z sign so that Z.real > 0
        #Z *= np.sign(Z.real)
        if sum(np.clip(Z.real,-10., 10.))<0: 
            Z *= -1
            #Z, Z_uo             = rt2z(s11, s12, init_sign=Z_init_sign)
else:
    N = np.zeros_like(freq)         # TODO why?
    Z = np.zeros_like(freq)
#}}}

## Detect resonances
losses = 1-abs(s11)**2-abs(s12)**2
loss_maxima = np.array(find_maxima(freq,losses))
#print "Detected loss maxima at frequencies:", loss_maxima
np.savetxt("last_found_modes.dat", loss_maxima)

## Get epsilon and mu
eps, mu = nz2epsmu(N, Z)

## Verify the results by back-calculating s11 and s12 to compare with the original values
s11backcalc, s12backcalc = nz2rt(freq, N, Z, d)

## --- Plotting to cartesian graphs -------------------------------------------- #{{{
plt.figure(figsize=(15,15))
xticks = np.arange(plot_freq_min, plot_freq_max, reasonable_ticks((plot_freq_max-plot_freq_min)/10))
xnumbers = [("%.2f"%(f/frequnit) if abs(f%reasonable_ticks((plot_freq_max-plot_freq_min)/10))<(frequnit/1000) else "") for f in xticks]
marker = "s" if (len(freq) < 20) else ""  # Use point markers for short data files
subplot_number = 4

## Plot reflection and transmission amplitudes
plt.subplot(subplot_number, 1, 1)
plt.plot(freq, s11amp, marker=marker, color="#AA4A00", label=u'$|s_{11}|$')
plt.plot(freq, s12amp, marker=marker, color="#004AAA", label=u'$|s_{12}|$')
plt.plot(freq, s12amp*1000, marker=marker, color="#00AA4A", label=u'$|s_{12}|*1000$')
plt.plot(freq, s12amp*100, marker=marker, color="#4AAA00", label=u'$|s_{12}|*100$')
plt.plot(freq, losses, color="#AAAAAA", label=u'loss')
if plot_expe and os.path.exists('r.dat'):
    rf, ry = np.loadtxt('r.dat', usecols=list(range(2)), unpack=True)
    plt.plot(rf, ry, lw=.5, ms=3, color='#AA8A40', marker='o') 
if plot_expe and os.path.exists('t.dat'):
    tf, ty = np.loadtxt('t.dat', usecols=list(range(2)), unpack=True)
    plt.plot(tf, ty, lw=.5, ms=3, color='#408AAA', marker='o') 

# - temporary -
if plot_expe and os.path.exists('../t00kVcm_Comsol.dat'):           ## XXX
    tf, ty = np.loadtxt('../t00kVcm_Comsol.dat', usecols=list(range(2)), unpack=True)
    plt.plot(tf*frequnit, ty, lw=2, color='#4A00AA', marker='o', alpha=.3, label='$|t_{0kV/cm}^{(Coms)}|$') 
if plot_expe and os.path.exists('../t90kVcm_Comsol.dat'):
    tf, ty = np.loadtxt('../t90kVcm_Comsol.dat', usecols=list(range(2)), unpack=True)
    plt.plot(tf*frequnit, ty, lw=2, color='#00AA4A', marker='s', alpha=.3, label='$|t_{90kV/cm}^{(Coms)}|$') 

## Verification of calculated data by calculating reflection and transmission again
plt.subplot(subplot_number, 1, 1) 
plt.plot(freq, abs(s11backcalc), color="#FA9962", label=u'$|s_{11FD}|$', ls='--')
plt.plot(freq, abs(s12backcalc), color="#6299FA", label=u'$|s_{12FD}|$', ls='--')
#plt.xticks(xticks, xnumbers); plt.minorticks_on(); plt.grid(1)

plt.ylabel(u"Amplitude"); plt.ylim((-0.1,1.1)); plt.xlim((plot_freq_min, plot_freq_max)) # XXX
plt.xticks(xticks, xnumbers); plt.minorticks_on();  plt.grid(True)
if legend_enable: plt.legend(loc="upper right"); 


#print '*********************************************'
#print 'R', s11amp[10]
#print 'T', s12amp[10]
#print 'L', 1-s11amp[10]**2-s12amp[10]**2
#print '*********************************************'

#for lm in loss_maxima: plt.axvspan(lm,lm+1e8, color='r')


## Plot r and t phase
# (Note: phase decreases with frequency, because meep uses the E=E0*exp(-i omega t) convention )
plt.subplot(subplot_number, 1, 2)
plt.plot(freq, np.unwrap(np.angle(s11))/pi, marker=marker, color="#AA4A00", label=u'$\\phi(s_{11})/\\pi$')
plt.plot(freq, np.unwrap(np.angle(s12))/pi, marker=marker, color="#004AAA", label=u'$\\phi(s_{12})/\\pi$')
#
#plt.plot(freq, np.unwrap(np.angle(s12))/pi + np.unwrap(np.angle(s11))/pi, marker=marker, color="#888AAA", label=u'$(\\phi(s_{11})+\\phi(s_{11}))/\\pi$')
#plt.plot(freq, np.unwrap(np.angle(s12))/pi - np.unwrap(np.angle(s11))/pi, marker=marker, color="#AA8A88", label=u'$(\\phi(s_{11})-\\phi(s_{11}))/\\pi$')
#plt.plot(freq, 2*np.unwrap(np.angle(s12))/pi + np.unwrap(np.angle(s11))/pi, marker=marker, color="#8A88AA", label=u'$(2\\phi(s_{11})+\\phi(s_{11}))/\\pi$')
#plt.plot(freq, 2*np.unwrap(np.angle(s12))/pi - np.unwrap(np.angle(s11))/pi, marker=marker, color="#8AAA88", label=u'$(2\\phi(s_{11})-\\phi(s_{11}))/\\pi$')
#plt.plot(freq, np.unwrap(np.angle(s12))/pi + 2*np.unwrap(np.angle(s11))/pi, marker=marker, color="#88AA8A", label=u'$(\\phi(s_{11})+2\\phi(s_{11}))/\\pi$')
#plt.plot(freq, np.unwrap(np.angle(s12))/pi - 2*np.unwrap(np.angle(s11))/pi, marker=marker, color="#AA888A", label=u'$(\\phi(s_{11})-2\\phi(s_{11}))/\\pi$')

# Optional: debugging curves(branch, sign, arg, anr, anl)
if len(freq)>2:
    #plt.plot(freq, N_debug[0]*.95, color="#dd0000", label=u"$br$", lw=1.6) 
    #plt.plot(freq, N_debug[1]*.90, color="#dd8800", label=u"$si$", lw=1.6) 
    #plt.plot(freq, N_debug[2].real, color="#00dd00", label=u"$arg^'$", lw=.6, ls='-') 
    #plt.plot(freq, N_debug[2].imag, color="#00dd00", label=u"$arg^{''}$", lw=.6, ls='--') 
    #plt.plot(freq, np.sign(N_debug[2].imag), color="#008800", label=u"sign$arg^{''}$", lw=.3, ls='-') 

    #plt.plot(freq, np.arccos(N_debug[2]).real, color="#0000dd", label=u"arccos$arg^'$", lw=1.6, ls='-') 
    #plt.plot(freq, np.log10(pi-np.arccos(N_debug[2]).real), color="#0000dd", label=u"arccos$arg^'$", lw=.6, ls='-') 
    #plt.plot(freq, np.arccos(N_debug[2]).imag, color="#0000dd", label=u"arccos$arg^{''}$", lw=1.6, ls='--') 
    #plt.plot(freq, np.log10(abs(N_debug[2].imag)), color="#000000", label=u"log$arg^{''}$", lw=.6, ls='--') 
    #plt.plot(freq, abs(N_debug[2] - (1+0j)), color="#0088dd", label=u"$|arg-1|$", lw=2, ls='-') 
    #plt.plot(freq, abs(N_debug[2] + (1+0j)), color="#8800dd", label=u"$|arg+1|$", lw=2, ls='-') 
    #plt.plot(freq, np.log10(abs(N_debug[2] - (1+0j))), color="#0088dd", label=u"", lw=1, ls='-') 
    #plt.plot(freq, np.log10(abs(N_debug[2] + (1+0j))), color="#8800dd", label=u"", lw=1, ls='-') 
    #plt.plot(freq, np.sign(N_debug[2].imag), color="#00dd00", label=u"$sgn arg^{''}$", lw=.6, ls=':') 
    plt.plot(freq, -np.ones_like(freq), color="k", label=u"", lw=.3, ls='-') 
    plt.plot(freq, np.ones_like(freq), color="k", label=u"", lw=.3, ls='-') 

    if autobranch:
        # Detection of key points in the spectrum (PBG boundaries, branch skips etc.)
        def find_maxima(x, y, minimum_value=.1):
            """ 
            Returns the x points where 
            1) y has a local maximum (i. e. dx/dy goes negative) AND 
            2) where y is above minimum_value 
            """
            d = y[1:-1] - y[0:-2]   ## naïve first derivative
            maxima = x[1:][np.sign(d[0:-2])-np.sign(d[1:-1]) + np.sign(y[2:-2]-minimum_value)==3]
            return maxima 
        def find_maxima_indices(x, y, minimum_value=.1):
            """ 
            Returns the x points where 
            1) y has a local maximum (i. e. dx/dy goes negative) AND 
            2) where y is above minimum_value 
            """
            d = y[1:-1] - y[0:-2]   ## naïve first derivative
            maximai = np.arange(1,len(x), dtype=np.dtype(np.int16))[np.sign(d[0:-2])-np.sign(d[1:-1]) + np.sign(y[2:-2]-minimum_value)==3]
            return maximai
        argPmin = find_maxima_indices(freq, -abs(N_debug[2] - (1+0j)), minimum_value=-np.inf)
        argNmin = find_maxima_indices(freq, -abs(N_debug[2] + (1+0j)), minimum_value=-np.inf)

        ## (todo) check: maybe required, maybe not
        #argNmax = find_maxima_indices(freq,  abs(N_debug[2] + (1+0j)), minimum_value=-np.inf)
        #plt.plot(freq[argNmax], np.zeros_like(argNmax), marker='o', color="#dd0000")
        #allindices = np.hstack([np.array([0]), argPmin, argNmin, argNmax])

        ## Concatenate &  sort all indices of interesting points
        allindices = np.hstack([np.array([0]), argPmin, argNmin])
        allindices.sort()
        ## Remove duplicate indices

        allindices = np.hstack([allindices[0], [x[0] for x in zip(allindices[1:],allindices[:-1]) if x[0]!=x[1]]])
        plt.plot(freq[allindices], np.zeros_like(allindices), marker='x', color="k")

        ## Scan through all photonic bands/bandgaps, seleting the correct N branch
        #print 'allindices', allindices
        #N_init_branch = 0
        #print 'N_init_sign', N_init_sign
        #N_init_sign = -1 
        #pN_uo = [0,0,0,0]
        pN_uo = [2*pi,2*pi,2*pi,0]
        det_branch = 0
        #for i in [0]:                           ## whole spectrum
                #i1 = 0
                #i2 = len(freq)-1
        for i in range(len(allindices)-1):     ## spectrum by chunks
            for q in (0,1):
                if q==0:
                    print 'LONG ',
                    i1 = allindices[i]
                    i2 = allindices[i+1]-1
                    #i2 = allindices[i+1]+1     ## .. works for 'q in [0]'
                else:
                    print 'SHORT',
                    i1 = allindices[i+1]-1
                    i2 = allindices[i+1]+1

                if i1>=i2: continue

                pfreq   = freq[i1:i2]
                if not q and pfreq[0] > 600e9: break
                pts = np.arange(10000)[i1:i2]; print pts[0], pts[-1],; print pfreq[0]/1e9,
                ps11    = s11[i1:i2]
                ps12    = s12[i1:i2]
                print 'start=', np.array(pN_uo)/pi,

                ## Plot oldschool N
                pN_uo_old = pN_uo
                pN, pN_uo, pN_debug    = rt2n(pfreq, ps11, ps12, d, init_branch=N_init_branch, init_sign=N_init_sign, uo=pN_uo)
                #if q!=0: pN_uo = pN_uo_old

                #print 'end=', np.array(pN_uo)/pi
                if i == 0:
                    try:
                        #print len(pN)
                        ii = 0
                        det_branch = np.round(2*np.real(pN[ii]*freq[ii]/c*d))
                        #print 'det_branch', det_branch
                    except:
                        pass
                    #print "N before correction", N[0:10]
                pN -= det_branch / (pfreq/c*d)/2
                plt.plot(pfreq, pN.real, lw=1.2, marker='o', markersize=2)
                #plt.plot(pfreq, pN.imag, lw=.8, ls='--')

                ## Plot oldschool UO
                #plt.plot(pfreq, np.ones_like(3pfreq)*pN_uo_old[0]/10, lw=3, c='#8888ff')
                #plt.plot(pfreq, np.ones_like(pfreq)*pN_uo_old[1]/10, lw=3, c='#88ff88', ls='-')
                #plt.plot(pfreq, np.ones_like(pfreq)*pN_uo_old[2]/10, lw=3, c='#ff8888', ls='-')
                #plt.plot(pfreq, np.ones_like(pfreq)*pN_uo_old[3]/10, lw=3, c='#88ffff', ls='-')




plt.ylabel(u"Phase"); None
plt.ylim((-15,15))
plt.xlim((plot_freq_min, plot_freq_max)) # XXX
#plt.xlim((00e9, 440e9))
plt.xticks(xticks, xnumbers); plt.minorticks_on(); plt.grid(True)
if legend_enable: plt.legend(); 


## Plot Z, N and figure-of-merit
plt.subplot(subplot_number, 1, 3)

if brillouin_boundaries:
    for i in range(1,4):
        plt.plot(freq, c/(2*freq*d)*i, color="#000000", label=u'', ls='-', lw=.5, alpha=.5)
        plt.plot(freq, -c/(2*freq*d)*i, color="#000000", label=u'', ls='-', lw=.5, alpha=.5)
if check_hilbert and len(freq)>1:
    import scipy.fftpack
    N[0] = N[1]  ## avoid NaN
    #np.kaiser(len(N), 5)
    N_KK = scipy.fftpack.hilbert(N.real + 1j*abs(N.imag)) / 1j  
    plt.plot(freq, np.real(N_KK), color="#FF9900", label=u"$N^{'}_{KK}$", alpha=1)
    plt.plot(freq, np.imag(N_KK), color="#FF9900", label=u'$N^{''}_{KK}$', ls='--', alpha=1)
    plt.plot(freq, np.real(N_KK)-np.real(N), color="#998800", label=u"$\\Delta N^{'}_{KK}$", alpha=.5, lw=3)
    plt.plot(freq, np.imag(N_KK)-np.imag(N), color="#998800", label=u'$\\Delta N^{''}_{KK}$', ls='--', alpha=.5, lw=3)

    Z[0] = Z[1]
    Z_KK = scipy.fftpack.hilbert(Z.real + 1j*Z.imag) / 1j   ## Why minus needed?
    #plt.plot(freq, np.real(Z_KK), color="#0099FF", label=u"$Z^{'}_{KK}$", alpha=.3)
    #plt.plot(freq, np.imag(Z_KK), color="#4499FF", label=u'$Z^{''}_{KK}$', ls='--', alpha=.3)
    DZr = np.real(Z_KK)-np.real(Z)
    DZi = np.imag(Z_KK)-np.imag(Z)
    #plt.plot(freq, DZr, color="#DDDD00", label=u"$\\Delta Z^{'}_{KK}$", alpha=.3)
    #plt.plot(freq, DZi, color="#DDDD44", label=u'$\\Delta Z^{''}_{KK}$', ls='--', alpha=.3)
    #plt.plot(freq[1:], (DZr[1:]+DZr[:-1])/2, color="#DDDD00", label=u"$\\Delta Z^{'}_{KK}$", alpha=.31)
    #plt.plot(freq[1:], (DZi[1:]+DZi[:-1])/2, color="#DDDD44", label=u'$\\Delta Z^{''}_{KK}$', ls='--', alpha=.31)


plt.plot(freq, np.real(N), color="#33AA00", label=u"$N$'")
plt.plot(freq, np.imag(N), color="#33AA33", label=u'$N$"', ls='--')

plt.plot(freq, np.real(Z), color="#0044DD", label=u"$Z$'")
plt.plot(freq, np.imag(Z), color="#4466DD", label=u'$Z$"', ls='--')

plt.plot(freq, np.log(-(np.real(N)/np.imag(N)))/np.log(10), 
    color="#FF9922", ls=":", label=u"$N$'$<0$ FOM")
plt.plot(freq, np.log((np.real(N)/np.imag(N)))/np.log(10), \
    color="#BB22FF", ls=":", label=u"$N$''$>0$ FOM")
plt.ylabel(u"Value"); 
plt.ylim((-5., 15.)); 
plt.xlim((plot_freq_min, plot_freq_max)); 
plt.xticks(xticks, xnumbers); plt.minorticks_on(); plt.grid(True)
if legend_enable: plt.legend(); 


## 4) Plot epsilon and mu
plt.subplot(subplot_number, 1, 4)

if find_plasma_frequency:
    try:
        from scipy.optimize import fsolve
        x, y = freq, eps.real
        estimates = x[np.where(np.diff(np.sign(y)))[0]]
        print "Plasma frequency (eps=0) at:", fsolve(lambda x0: np.interp(x0, x, y),  estimates)
    except:
        print "Plasma frequency (epsilon(f) == 0) detection failed"
     
plt.xlabel(u"Frequency [%s]" % frequnitname) 
if plot_expe and os.path.exists('eps.dat'):
    tf, ty = np.loadtxt('eps.dat', usecols=list(range(2)), unpack=True)
    plt.plot(tf*frequnit, ty, lw=0, color='#AA0088', marker='o') ## XXX
    plt.plot(tf*frequnit, -ty, lw=0, color='#AA8888', marker='s') ## XXX
    #plt.plot(tf     , ty, lw=0, color='#AA0088', marker='o') ## XXX
if plot_expe and os.path.exists('mu.dat'):
    tf, ty = np.loadtxt('mu.dat', usecols=list(range(2)), unpack=True)
    plt.plot(tf*frequnit, ty, lw=0, color='#AA8800', marker='o') ## XXX
    plt.plot(tf*frequnit, -ty, lw=0, color='#AA8888', marker='s') ## XXX
    #plt.plot(tf     , ty, lw=0, color='#AA0088', marker='o') ## XXX

if check_hilbert and len(freq)>1:
    import scipy.fftpack
    eps[0] = 0  ## avoid NaN
    eps_KK = scipy.fftpack.hilbert(eps.real + 1j*abs(eps.imag)) / 1j  
    plt.plot(freq, np.real(eps_KK), color="#FF9900", label=u"$eps^{'}_{KK}$", alpha=.5)
    plt.plot(freq, np.imag(eps_KK), color="#FF9900", label=u"$eps^{''}_{KK}$", ls='--', alpha=.5)
    plt.plot(freq, np.real(eps_KK)-np.real(eps), color="#FF0099", label=u"$eps^{'}_{KK}$", alpha=.5)
    plt.plot(freq, np.imag(eps_KK)-np.imag(eps), color="#FF0099", label=u"$eps^{''}_{KK}$", ls='--', alpha=.5)

    mu[0] = 0
    mu_KK = scipy.fftpack.hilbert(N.real + 1j*abs(N.imag)) / 1j  
    plt.plot(freq, np.real(mu_KK), color="#0099FF", label=u"$mu^{'}_{KK}$", alpha=.5)
    plt.plot(freq, np.imag(mu_KK), color="#4499FF", label=u'$mu^{''}_{KK}$', ls='--', alpha=.5)
    plt.plot(freq, np.real(mu_KK)-np.real(mu), color="#0099FF", label=u"$mu^{'}_{KK}$", alpha=.5)
    plt.plot(freq, np.imag(mu_KK)-np.imag(mu), color="#4499FF", label=u'$mu^{''}_{KK}$', ls='--', alpha=.5)
plt.plot(freq, np.real(eps), color="#AA0088", label=u"$\\varepsilon_{eff}$'")
plt.plot(freq, np.imag(eps), color="#FF22DD", label=u'$\\varepsilon_{eff}$"', ls='--')
plt.plot(freq, np.real(mu),  color="#BB8800", label=u"$\\mu_{eff}$'")
plt.plot(freq, np.imag(mu),  color="#DDAA00", label=u'$\\mu_{eff}$"', ls='--')
plt.ylabel(u"Value"); plt.ylim((-10.,10.)); 
plt.yscale('symlog', linthreshy=10.); 
plt.xlim((plot_freq_min, plot_freq_max))
plt.xticks(xticks, xnumbers); plt.minorticks_on(); plt.grid(True)
if legend_enable: plt.legend(); 

## Final plotting 
plt.savefig(last_simulation_name+".png", bbox_inches='tight')
#}}}
## --- Plotting to dispersion curves (k-omega) -------------------------------------------- #{{{
if plot_bands and not os.path.exists("band"): os.mkdir("band")
if plot_bands and os.path.isdir("band"):
    plt.figure(figsize=(8,8))
    plt.plot(np.arcsin(np.sin(np.real(N*freq*d/c) * pi)) / pi, freq, color="#33AA00", label=u"$k$'")
    plt.plot(np.imag(N*freq*d/c), freq, color="#33AA33", label=u'$\\kappa$', ls='--')

    ## Detection of bandgap: ratio of the real to the imaginary part of complex wavenumber
    ## we will use the sin() of the real part of k so that it does not matter in which Brillouin zone it is
    try:
        realpart = np.arcsin(np.sin(pi * 2*np.real(N*freq/c*d)))
        imagpart = np.abs(np.imag(N*freq/c*d))
        pbg_indicator = np.sign(abs(realpart) - abs(imagpart))
        ## starts and ends of band-gap
        pbg_starts = np.interp(np.where(pbg_indicator[1:] < pbg_indicator[0:-1]), range(len(freq)), freq)[0]
        pbg_ends   = np.interp(np.where(pbg_indicator[1:] > pbg_indicator[0:-1]), range(len(freq)), freq)[0] 
        ## Fix the un-started and un-ended bandgaps
        if len(pbg_starts) < len(pbg_ends): pbg_starts = np.concatenate([np.array([0]), pbg_starts])
        if len(pbg_starts) > len(pbg_ends): pbg_starts = pbg_starts[:-1]
        for start, end in np.vstack([pbg_starts, pbg_ends]).T:
            plt.axhspan(start, end, color='#FFDD00', alpha=.1)
    except:
        print "Bandgap detection failed"

    plt.ylabel(u"frequency"); 
    plt.xlabel(u"wavenumber $ka/\\pi$"); 
    plt.xlim((-.5, .5)); 
    plt.xticks(xticks, xnumbers); plt.minorticks_on(); 
    plt.grid(True)
    if legend_enable: plt.legend(loc="upper right"); 

    ## Final plotting 
    splitpath = os.path.split(last_simulation_name)
    outfile = os.path.join(splitpath[0], "band", splitpath[1]+"_band.png")
    plt.savefig(outfile, bbox_inches='tight')
#}}}
## --- Nice plotting to PDF ----------------------------------------------------------------------------------#{{{
#if plot_publi and not os.path.exists("publi"): os.mkdir("publi") TODO remove
if plot_publi:
    if not os.path.exists("publi"): os.mkdir("publi")
    #matplotlib.rc('text', usetex=True)
    #matplotlib.rc('text.latex', preamble = \
            #'\usepackage{amsmath}, \usepackage{yfonts}, \usepackage{txfonts}, \usepackage{lmodern},')

    matplotlib.rc('text', usetex=True)
    matplotlib.rc('font', size=14)
    matplotlib.rc('text.latex', preamble = \
            '\usepackage{amsmath}, \usepackage{palatino},\usepackage{upgreek}')
    matplotlib.rc('font',**{'family':'serif','serif':['palatino, times']})  ## select fonts

    fig = plt.figure(figsize=(12,10))
    fig.subplots_adjust(left=.05, bottom=.05, right=.99, top=.99, wspace=.0, hspace=.0) ## XXX

    publi_toplot = {'rt':1, 'N':1, 'eps':1, 'mu':1, 'Z':0} ## select parameters by setting 1 or 0
    publi_plot_Brillouin = True
    publi_use_grid = False

    subplot_index  = 1
    subplot_count  = sum(publi_toplot.values())
    ## ---- r, t -----
    if publi_toplot['rt']:
        ax= plt.subplot(subplot_count, 1, subplot_index)
        #plt.title(u"Dielectric spheres $r=%d\\;\\upmu$m" % 25) 
        #plt.title(u"Dielectric spheres in wire mesh") 
        plt.title(u"Wire mesh") 
        ax.label_outer()
        plt.grid(publi_use_grid)
        plt.plot(freq, s11amp, marker=marker, color="#880000", label=u'$|r|$', lw=1)
        plt.plot(freq, s12amp, marker=marker, color="#0088ff", label=u'$|t|$', lw=1)
        plt.ylabel(u"Amplitude"); 
        if plot_expe and os.path.exists('t.dat'):
            tf, ty = np.loadtxt('t.dat', usecols=list(range(2)), unpack=True)
            plt.plot(tf*frequnit, ty, lw=0, color='#004AAA', marker='o', ms=2, label=u'$|t|$ exp') 
        subplot_index += 1
        plt.xticks(xticks, xnumbers); plt.minorticks_on(); 
        plt.xlim((plot_freq_min, plot_freq_max)); plt.ylim((0,1.)); plt.legend(loc='lower left');

    ## Todo allow plotting phase! (And in the 'cartesian' plot, too)

    ## ---- N -----
    if publi_toplot['N']:
        ax = plt.subplot(subplot_count, 1, subplot_index)
        ax.label_outer()
        plt.grid(publi_use_grid)
        plt.ylabel(u"Index of refraction  $N_{\\text{eff}}$"); 

        if publi_plot_Brillouin:
            for ii in np.arange(-10, 10):
                plt.plot(freq, ii*c/freq/d, color="#000000", label=u"", lw=.2)
                plt.plot(freq, (ii+.5)*c/freq/d, color="#777777", label=u"", lw=.2)

        #TODO if plot_expe and os.path.exists('k.dat'):
            #tf, ty = np.loadtxt('t.dat', usecols=list(range(2)), unpack=True)
            #plt.plot(tf*frequnit, ty, lw=0, color='#004AAA', marker='o', ms=2, label=u'$|t|$ exp') 

        plt.plot(freq, np.real(N), color="#448800", label=u"$N'$")
        plt.plot(freq, np.imag(N), color="#448800", label=u"$N''$", ls='--')
        if check_hilbert and len(freq)>1:
            plt.plot(freq, np.real(N_KK), color="#dd88aa", label=u"")
            plt.plot(freq, np.imag(N_KK), color="#dd88aa", label=u"", ls='--')
        plt.xticks(xticks, xnumbers); plt.minorticks_on()
        plt.xlim((plot_freq_min, plot_freq_max)); plt.ylim((-5,5)); plt.legend(loc='lower left'); 
        subplot_index += 1

    ## ----- EPS -----
    if publi_toplot['eps']:
        ax = plt.subplot(subplot_count, 1, subplot_index)
        ax.label_outer()
        plt.grid(publi_use_grid)
        plt.ylabel(u"Permittivity $\\varepsilon_{\\text{eff}}$") 
        plt.plot(freq, np.real(eps), color="#660044", label=u"$\\varepsilon'$")
        plt.plot(freq, np.imag(eps),       color="#660044", label=u"$\\varepsilon''$", ls='--')

        ## optional: Drude model
        #plt.plot(freq, 1-(1100e9/freq)**2,       color="#888888", label=u"$1-\\frac{f_p^2}{f^2}$", ls='-') 

        plt.xticks(xticks, xnumbers); plt.minorticks_on()
        plt.xlim((plot_freq_min, plot_freq_max)); plt.ylim((-5.,5.)); plt.legend(loc='lower left'); 
        subplot_index += 1

    ## ----- MU -----
    if publi_toplot['mu']:
        ax = plt.subplot(subplot_count, 1, subplot_index)
        ax.label_outer()
        plt.grid(publi_use_grid)
        plt.ylabel(u"Permeability $\\mu_{\\text{eff}}$"); 
        plt.plot(freq, np.real(mu), color="#663300", label=u"$\\mu'$")
        plt.plot(freq, np.imag(mu), color="#663300", label=u"$\\mu''$", ls='--')
        plt.xticks(xticks, xnumbers); plt.minorticks_on(); 
        plt.xlim((plot_freq_min, plot_freq_max)); 
        plt.ylim((-5.,5.)); 
        plt.legend(loc='lower left'); 
        subplot_index += 1

    ### ----- Z -----
    if publi_toplot['Z']:
        ax = plt.subplot(subplot_number, 1, subplot_index)
        ax.label_outer()
        plt.ylabel(u"Impedance"); plt.ylim((-2.,4.))
        plt.plot(freq, np.real(Z), color="#004488", label=u"$Z'$")
        plt.plot(freq, np.imag(Z), color="#004488", label=u"$Z''$", ls='--')
        plt.xticks(xticks, xnumbers); plt.minorticks_on(); 
        plt.xlim((plot_freq_min, plot_freq_max));  plt.legend(loc=(.03,.6));
        subplot_index += 1

    plt.xlabel(u"Frequency [%s]" % frequnitname) 
    splitpath = os.path.split(last_simulation_name)
    outfile = os.path.join(splitpath[0], "publi", splitpath[1]+"_publi.pdf")
    plt.savefig(outfile, bbox_inches='tight')
    #}}}

## --- Save data to effparam.dat (to ./effparam/*dat) ------------------------------------------#{{{
if savedat:
    splitpath = os.path.split(last_simulation_name)
    if not os.path.exists("effparam"): os.mkdir("effparam") ## FIXME - fails if processing file outside pwd?
    savedatfile = os.path.join(splitpath[0], "effparam", splitpath[1]+"_effparam.dat")

    ## Copy parameters - load header 
    header = ""
    with open(last_simulation_name+".dat") as datafile:
        for line in datafile:
            if (line[:1]=="#") and (not "olumn" in line): header+=line
    with open(savedatfile, "w") as outfile:
        ## Copy parameters - write header 
        outfile.write(header)
        ## Write column headers
        outfile.write("#x-column freq\n#column |r|\n#column r phase\n#column |t|\n#column t phase\n" + \
                    "#column real N\n#column imag N\n#column real Z\n#column imag Z\n" + \
                    "#column real eps\n#column imag eps\n#column real mu\n#column imag mu\n")
        ## Write column data
        np.savetxt(outfile, zip(freq, s11amp, s11phase, s12amp, s12phase, 
                    N.real, N.imag, Z.real, Z.imag, eps.real, eps.imag, mu.real, mu.imag), fmt="%.8e")
#}}}
## --- Plot polar ------------------------------------------------------------#{{{
if plot_polar and not os.path.exists("polar"): os.mkdir("polar")
if plot_polar and os.path.isdir("polar"):
    ## Truncate the arrays (optional)
    #(d0,d1) = np.interp((500e9, 650e9), freq, range(len(freq)))
    #(freq, s11, s12, N, Z, eps, mu, arg, argLog) = \
            #map(lambda a: a[int(d0):int(d1)], (freq, s11, s12, N, Z, eps, mu, arg, argLog))

    print "Plotting polar..."
    from matplotlib.collections import LineCollection
    lims={"s11":(-1,1), "s12":(-1,1), "N":(-10,10), "Z":(-5,5), 
            "mu":(-10,10), "eps":(-10,10), "arg":(-3,3), "argLog":(-10,10) }
    datalist=(s11, s12, N, Z, eps, mu, arg, argLog)

    plotlabels=("s11", "s12", "N", "Z", "eps", "mu", "arg", "argLog")

    freqlabels = np.append(loss_maxima[loss_maxima<plot_freq_max], freq[-1])

    fig = plt.figure(figsize=(11,22))
    subplot_number = len(datalist)
    for (subpl, data, plotlabel) in zip(range(subplot_number), datalist, plotlabels):
        plt.subplot(4,2,subpl+1)
        if plotlabel.startswith('s'):
            plt.plot(np.sin(np.linspace(0,2*pi)), np.cos(np.linspace(0,2*pi)), c='#888888')
            plt.plot(np.sin(np.linspace(0,2*pi))/2+.5, np.cos(np.linspace(0,2*pi))/2, c='#aaaaaa')
            plt.plot(np.sin(np.linspace(0,2*pi))+1, np.cos(np.linspace(0,2*pi))+1, c='#aaaaaa')
            plt.plot(np.sin(np.linspace(0,2*pi))+1, np.cos(np.linspace(0,2*pi))-1, c='#aaaaaa')

        x = data.real; y = data.imag
        t = np.linspace(0, 10, len(freq))
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        lc = LineCollection(segments, cmap=plt.get_cmap('jet'), norm=plt.Normalize(0, 10))
        lc.set_array(t)
        lc.set_linewidth(2)

        plt.gca().add_collection(lc)

        ## Add black points to every xtick
        xpoints = np.interp(xticks, freq, x)
        ypoints = np.interp(xticks, freq, y)
        for xpoint, ypoint in zip(xpoints, ypoints):
            plt.plot(xpoint, ypoint, marker="o", markersize=3, color="#000000", label='')

        ## Annotate resonant frequencies
        xpoints = np.interp(freqlabels, freq, x.real)
        ypoints = np.interp(freqlabels, freq, y.real)
        freqlabelstxt = [("%d" % (fr*1000/frequnit)) for fr in freqlabels]
        for label, xpoint, ypoint in zip(freqlabelstxt, xpoints, ypoints):
            plt.annotate(label, xy = (xpoint, ypoint), xytext = (-10, 10),
                textcoords = 'offset points', ha = 'right', va = 'bottom',
                bbox = dict(boxstyle = 'round,pad=.15', fc = 'white', alpha = 0.5),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            plt.plot(xpoint, ypoint, marker="o", markersize=2, color="#000000", label='')

        lim = lims[plotlabel]
        plt.xlim(lim); plt.ylim(lim); plt.grid(True); plt.title(plotlabel)
    ## Final plotting
    splitpath = os.path.split(last_simulation_name)
    outfile = os.path.join(splitpath[0], "polar", splitpath[1]+"_polar.png")
    plt.savefig(outfile, bbox_inches='tight')

#}}}


