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
import sys, os, re, matplotlib, argparse
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, fmin
from scipy.constants import pi, c

## == User settings for postprocessing and plotting == 
parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
## Effective parameter retrieval options
parser.add_argument('--N_init_sign',                    type=int,   default=1)
parser.add_argument('--N_init_branch',                  type=int,   default=0)
parser.add_argument('--Z_init_sign',                    type=int,   default=-1)
parser.add_argument('--padding',                        type=float, default=0., help='')
parser.add_argument('--numstabcoef',                    type=float, default=.9997, help='slight artificial "absorption" to keep algorithm numerically stable')
parser.add_argument('--autocorrect_signs',              type=int,   default=1, help='shall enforce the positive value of imag N?')
parser.add_argument('--autobranch_sampler_position',    type=float, default=0.03)
parser.add_argument('--autobranch',                     type=int,   default=0, help='automatepsic selection of the branch')
parser.add_argument('--autocorrect_signs_pointwise',    type=int,   default=0, help='shall finally enforce the positive value of imag N and real Z in every point?')
parser.add_argument('--savedat',                        type=int,   default=1, help='created directory "effparam" and saves all params to an ascii file with header')
## Postprocessing
parser.add_argument('--find_plasma_frequency',          type=int,   default=0, help='find frequencies where epsilon crosses zero')
parser.add_argument('--detect_loss_maxima',             type=int,   default=0, help='')
## Plotting mode choice
parser.add_argument('--plot_polar',                     type=int,   default=0, help='plots results to polar graphs for diagnostics')
parser.add_argument('--plot_bands',                     type=int,   default=0, help='plots index of refraction as dispersion curves (k-omega)')
## Plotting tuning
parser.add_argument('--frequnit',                       type=int,   default=1e12)
parser.add_argument('--frequnitname',                   type=str,   default='THz')
parser.add_argument('--check_hilbert',                  type=int,   default=0, help='Verifies if Kramers-Kronig relations hold for N')  ###XXX?
parser.add_argument('--legend_enable',                  type=int,   default=1, help='show legend in the plot')
parser.add_argument('--brillouin_boundaries',           type=int,   default=1, help='Plots thin lines where the N would exceed the allowed range for 0-th Bloch mode')
parser.add_argument('--plot_expe',                      type=int,   default=1, help='if "r.dat", "t.dat", "N.dat", "Z.dat", "eps.dat" or "mu.dat" available, overlay them')
parser.add_argument('--plot_freq_min',                  type=float, default=np.nan, help='')
parser.add_argument('--plot_freq_max',                  type=float, default=np.nan, help='if None, decide from the input file header')
parser.add_argument('--plot_weak_transmission',         type=int,   default=0)
parser.add_argument('filenames',    type=str,   nargs='?', help='DAT files to be processed')
args = parser.parse_args()

np.seterr(all='ignore')      ## do not print warnings for negative-number logarithms etc.


## == Auxiliary functions ==
def get_simulation_name(argindex=1): #{{{
    """Get the name of the last simulation run.

    Priority: 1) parameter, 2) last_simulation_name.dat, 3) working directory"""
    cwd = os.getcwd()
    print args.filenames
    if args.filenames  and __name__ == "__main__": 
        print "Parameter passed:", args.filenames
        last_simulation_name = args.filenames
    elif os.path.exists(os.path.join(cwd, 'last_simulation_name.dat')):
        print "Loading from", os.path.join(cwd, 'last_simulation_name.dat')
        last_simulation_name = os.path.join(cwd, open(os.path.join(cwd, 'last_simulation_name.dat'),'r').read().strip())
        print 'last_simulation_name',last_simulation_name
    else:
        print "Error: No input file provided and 'last_simulation_name.dat' not found!"
        last_simulation_name = cwd
    if (last_simulation_name[-4:] == ".dat"): last_simulation_name = last_simulation_name[:-4] # strip the .dat extension
    return  last_simulation_name
#}}}
def load_rt(filename): #{{{
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
    #with open(filename+'.dat') as datafile:
        #for line in datafile:
            #if line[0:1] in "0123456789": break         # end of file header
            #value = line.replace(",", " ").split()[-1]  # the value of the parameter will be separated by space or comma
            #if ("cellsize" in line) and (cellsize == None): cellsize = float(value)
            #if ("cellnumber" in line): cellnumber = float(value)
            #if ("plot_freq_min" in line) and (plot_freq_min == 0): plot_freq_min = float(value)
            #if ("plot_freq_max" in line) and (plot_freq_max == np.infty): plot_freq_max = float(value)
            #if ("param padding" in line) and (padding == None): padding = float(value)

    ## Load data columns
    (freq, s11amp, s11phase, s12amp, s12phase) = \
            map(lambda a: np.array(a, ndmin=1), np.loadtxt(filename+".dat", unpack=True)) 
    
    ## Limit the frequency range to what will be plotted (recommended) TODO wrong approach
    #XXX RM if truncate and len(freq)>1:
        #(d0,d1) = np.interp((plot_freq_min, plot_freq_max), freq, range(len(freq)))
        #(freq, s11amp, s11phase, s12amp, s12phase) = \
                #map(lambda a: a[int(d0):int(d1)], (freq, s11amp, s11phase, s12amp, s12phase))
    return freq, s11amp, s11phase, s12amp, s12phase
#}}}
def get_param(filename):             ## Load header to the 'parameters' dictionary#{{{
    parameters = {}
    with open(filename) as datafile:
        for line in datafile:
            if (line[0:1] in '0123456789') or ('column' in line.lower()): break    # end of parameter list
            ## key-value separator is either ',' or '='; take the left word from it as the param name, and everything on the right as the param value
            left, value = line.replace(',', '=').strip().split('=', 1)
            key = left.split()[-1]
            try: value = float(value) ## Try to convert to float, if possible
            except: pass                ## otherwise keep as string
            parameters[key] = value
    return parameters
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
    1) This optimum shift shall hold for all simulations with one or more unit cellnumber.
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

## == Homogenisation functions (efficient processing whole np.array at once) ==
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

          Putting the layers far apart alleviates this for 2 cellnumber: can it be related to higher-order Bloch modes? 

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

## --- Preparation of data ------------------------------------ # {{{
## Get reflection and transmission data, prepare its parameters
last_simulation_name = get_simulation_name()
print 'last_simulation_name', last_simulation_name
freq, s11amp, s11phase, s12amp, s12phase = load_rt(last_simulation_name)
params  = get_param(last_simulation_name+'.dat')

if 'cellsize' in params.keys():   cellsize = params['cellsize']
else:                             print "Warning, `cellsize' parameter not specified, effparam retrieval wrong"; cellsize = 1 

if 'cellnumber' in params.keys(): cellnumber = params['cellnumber']
else:                             print "Warning, `cellnumber' parameter not specified, defaulting to 1"; cellnumber = 1 

d = cellsize * cellnumber         ## total thickness of the structure

if 'padding' in params.keys():    padding = params['padding']
else:                             print "Warning, `padding' parameter not specified, defaulting to 0"; padding = 1 

if not args.plot_freq_min is np.nan: plot_freq_min = args.plot_freq_min
elif 'plot_freq_min' in params:      plot_freq_min = params['plot_freq_min']
else:                                plot_freq_min = np.min(freq)

if not args.plot_freq_max is np.nan: plot_freq_max = args.plot_freq_max
elif 'plot_freq_max' in params:      plot_freq_max = params['plot_freq_max']
else:                                plot_freq_max = np.max(freq)
# }}}
## --- Calculation of effective parameters -------------------- # {{{
## Convert to complex numbers and compensate for the additional padding of the monitor planes
s11 = shiftmp(freq, polar2complex(s11amp, s11phase), padding*np.ones_like(freq)) * args.numstabcoef
s12 = shiftmp(freq, polar2complex(s12amp, s12phase), padding*np.ones_like(freq)) * args.numstabcoef

## Build the debug plots
arg = (1+0j-s11**2+s12**2)/2/(s12)
argLog = np.e**(1j*np.angle(arg))*np.log(1+abs(arg)) ## radially shrinked graph to see the topology

## Calculate N, Z and try to correct the signs
if len(freq)>2:
    N, N_uo, N_debug    = rt2n(freq, s11, s12, d, init_branch=args.N_init_branch, init_sign=args.N_init_sign)
    #print "N before correctio1", N[0:10]
    Z, Z_uo             = rt2z(s11, s12, init_sign=args.Z_init_sign)
    if args.autocorrect_signs: 

        ## Fix N sign so that N.imag > 0 
        nimag = np.sum(np.clip(N.imag,-10., 10.))
        if nimag<0: N *= -1

        ## Take a sample
        ii = int(float(len(N)) * args.autobranch_sampler_position)
        sampleR = np.real(N[ii])
        sampleI = np.imag(N[ii])
        #print 'sampleR, sampleI =', sampleR, sampleI

        #if ii == -0: N *= -1  ## alternative way of fixing sign?

        ## Fix N branch so that N.real does not diverge  at low frequencies
        branch_selector = 2*sampleR*freq[ii]/c*d
        det_branch = np.round(branch_selector)
        if (branch_selector-det_branch) < 0:
            N *= -1
            det_branch *= -1
            branch_selector *= -1
        #print 'branch_selector, det_branch, diff', branch_selector, det_branch, (branch_selector-det_branch)
        N -= det_branch / (freq/c*d)/2

        ## Fixing Z sign so that Z.real > 0
        #if sum(np.clip(Z.real,-10., 10.))<0: 
        if np.real(Z[ii])<0: 
            Z *= -1

    if args.autocorrect_signs_pointwise: 
        N *= np.sign(N.imag)
        Z *= np.sign(Z.imag)
else:
    N = np.zeros_like(freq) 
    Z = np.zeros_like(freq)

## Detect resonances
losses = 1-abs(s11)**2-abs(s12)**2
loss_maxima = np.array(find_maxima(freq,losses))
#print "Detected loss maxima at frequencies:", loss_maxima
np.savetxt("last_found_modes.dat", loss_maxima)

## Get epsilon and mu
eps, mu = nz2epsmu(N, Z)

## Verify the results by back-calculating s11 and s12 to compare with the original values
s11backcalc, s12backcalc = nz2rt(freq, N, Z, d)
# }}}

## --- Plotting to cartesian graphs -------------------------------------------- #{{{
plt.figure(figsize=(15,15))
marker = "s" if (len(freq) < 20) else ""  # Use point markers for short data files
subplot_number = 4

## Plot reflection and transmission amplitudes
plt.subplot(subplot_number, 1, 1)
plt.plot(freq, s11amp, marker=marker, color="#AA4A00", label=u'$|s_{11}|$')
plt.plot(freq, s12amp, marker=marker, color="#004AAA", label=u'$|s_{12}|$')
if args.plot_weak_transmission:
    plt.plot(freq, s12amp*1000, marker=marker, color="#00AA4A", label=u'$|s_{12}|*1000$')
    plt.plot(freq, s12amp*100, marker=marker, color="#4AAA00", label=u'$|s_{12}|*100$')
plt.plot(freq, losses, color="#AAAAAA", label=u'loss')
if args.plot_expe and os.path.exists('r.dat'):
    rf, ry = np.loadtxt('r.dat', usecols=list(range(2)), unpack=True)
    plt.plot(rf, ry, lw=.5, ms=3, color='#AA8A40', marker='o') 
if args.plot_expe and os.path.exists('t.dat'):
    tf, ty = np.loadtxt('t.dat', usecols=list(range(2)), unpack=True)
    plt.plot(tf, ty, lw=.5, ms=3, color='#408AAA', marker='o') 

## Verification of calculated data by calculating reflection and transmission again
plt.subplot(subplot_number, 1, 1) 
plt.plot(freq, abs(s11backcalc), color="#FA9962", label=u'$|s_{11FD}|$', ls='--')
plt.plot(freq, abs(s12backcalc), color="#6299FA", label=u'$|s_{12FD}|$', ls='--')

plt.ylabel(u"Amplitude"); plt.ylim((-0.1,1.1)); plt.xlim((plot_freq_min, plot_freq_max)) # XXX
if args.legend_enable: plt.legend(loc="upper right"); 

## Plot r and t phase
plt.subplot(subplot_number, 1, 2)
plt.plot(freq, np.unwrap(np.angle(s11))/pi, marker=marker, color="#AA4A00", label=u'$\\phi(s_{11})/\\pi$')
plt.plot(freq, np.unwrap(np.angle(s12))/pi, marker=marker, color="#004AAA", label=u'$\\phi(s_{12})/\\pi$')
plt.ylabel(u"Phase")
plt.ylim((-15,15))
plt.xlim((plot_freq_min, plot_freq_max)) # XXX
if args.legend_enable: plt.legend(); 


## Plot Z, N and figure-of-merit
plt.subplot(subplot_number, 1, 3)

if args.brillouin_boundaries:
    for i in range(1,4):
        plt.plot(freq, c/(2*freq*d)*i, color="#000000", label=u'', ls='-', lw=.5, alpha=.5)
        plt.plot(freq, -c/(2*freq*d)*i, color="#000000", label=u'', ls='-', lw=.5, alpha=.5)
if args.check_hilbert and len(freq)>1:
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

plt.plot(freq, np.real(N), color="#33AA00", label=u"$N$'")
plt.plot(freq, np.imag(N), color="#33AA33", label=u'$N$"', ls='--')

#plt.plot(freq, np.real(Z), color="#0044DD", label=u"$Z$'")
#plt.plot(freq, np.imag(Z), color="#4466DD", label=u'$Z$"', ls='--')

#plt.plot(freq, np.log(-(np.real(N)/np.imag(N)))/np.log(10), 
    #color="#FF9922", ls=":", label=u"$N$'$<0$ FOM")
#plt.plot(freq, np.log((np.real(N)/np.imag(N)))/np.log(10), \
    #color="#BB22FF", ls=":", label=u"$N$''$>0$ FOM")
plt.ylabel(u"Value"); 
plt.ylim((-3., 6.)); 
plt.xlim((plot_freq_min, plot_freq_max)); 
plt.grid(True)
if args.legend_enable: plt.legend(); 


## 4) Plot epsilon and mu
plt.subplot(subplot_number, 1, 4)

if args.find_plasma_frequency:
    try:
        from scipy.optimize import fsolve
        x, y = freq, eps.real
        estimates = x[np.where(np.diff(np.sign(y)))[0]]
        print "Plasma frequency (eps=0) at:", fsolve(lambda x0: np.interp(x0, x, y),  estimates)
    except:
        print "Plasma frequency (epsilon(f) == 0) detection failed"
     
plt.xlabel(u"Frequency [%s]" % args.frequnitname) 
if args.plot_expe and os.path.exists('eps.dat'):
    tf, ty = np.loadtxt('eps.dat', usecols=list(range(2)), unpack=True)
    plt.plot(tf*args.frequnit, ty, lw=0, color='#AA0088', marker='o') ## XXX
    plt.plot(tf*args.frequnit, -ty, lw=0, color='#AA8888', marker='s') ## XXX
    #plt.plot(tf     , ty, lw=0, color='#AA0088', marker='o') ## XXX
if args.plot_expe and os.path.exists('mu.dat'):
    tf, ty = np.loadtxt('mu.dat', usecols=list(range(2)), unpack=True)
    plt.plot(tf*args.frequnit, ty, lw=0, color='#AA8800', marker='o') ## XXX
    plt.plot(tf*args.frequnit, -ty, lw=0, color='#AA8888', marker='s') ## XXX 
    #plt.plot(tf     , ty, lw=0, color='#AA0088', marker='o') ## XXX 
if args.check_hilbert and len(freq)>1:
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
#plt.yscale('symlog', linthreshy=.1); 
plt.xlim((plot_freq_min, plot_freq_max))
plt.grid(True)
if args.legend_enable: plt.legend(); 

## Final plotting 
plt.savefig(last_simulation_name+".png", bbox_inches='tight')
#}}}
## --- Plotting to dispersion curves (k-omega) -------------------------------------------- #{{{
if args.plot_bands and not os.path.exists("band"): os.mkdir("band")
if args.plot_bands and os.path.isdir("band"):
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
    plt.grid(True)
    if args.legend_enable: plt.legend(loc="upper right"); 

    ## Final plotting 
    splitpath = os.path.split(last_simulation_name)
    outfile = os.path.join(splitpath[0], "band", splitpath[1]+"_band.png")
    plt.savefig(outfile, bbox_inches='tight')
#}}}

## --- Save data to effparam.dat (to ./effparam/*dat) ------------------------------------------#{{{
if args.savedat:
    splitpath = os.path.split(last_simulation_name)
    if not os.path.exists("effparam"): os.mkdir("effparam") ## FIXME - fails if processing file outside pwd?
    args.savedatfile = os.path.join(splitpath[0], "effparam", splitpath[1]+"_effparam.dat")

    ## Copy parameters - load header 
    header = ""
    with open(last_simulation_name+".dat") as datafile:
        for line in datafile:
            if (line[:1]=="#") and (not "olumn" in line): header+=line
    with open(args.savedatfile, "w") as outfile:
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
if args.plot_polar and not os.path.exists("polar"): os.mkdir("polar")
if args.plot_polar and os.path.isdir("polar"):
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

        ## Add black points to every xtick DEFUNCT
        #xpoints = np.interp(xticks, freq, x)
        #ypoints = np.interp(xticks, freq, y)
        #for xpoint, ypoint in zip(xpoints, ypoints):
            #plt.plot(xpoint, ypoint, marker="o", markersize=3, color="#000000", label='')

        ## Annotate resonant frequencies
        xpoints = np.interp(freqlabels, freq, x.real)
        ypoints = np.interp(freqlabels, freq, y.real)
        freqlabelstxt = [("%d" % (fr*1000/args.frequnit)) for fr in freqlabels]
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
    outfile = os.path.join(splitpath[0], "polar", splitpath[1]+"_polar.pdf")
    plt.savefig(outfile, bbox_inches='tight')

#}}}

