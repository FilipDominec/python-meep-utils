import time, sys, os
import numpy as np
from scipy.constants import c, epsilon_0, mu_0

import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep

maxftmp = 6e9

## Get the name of the last simulation run 
cwd = os.getcwd()
argindex=1
if len(sys.argv)>argindex and sys.argv[argindex] != "-"  and __name__ == "__main__": 
    print "Parameter passed:", sys.argv[argindex]
    last_simulation_name = sys.argv[argindex]
elif os.path.exists(os.path.join(cwd, 'last_simulation_name.dat')):
    print "Loading from", os.path.join(cwd, 'last_simulation_name.dat')
    last_simulation_name = os.path.join(cwd, open(os.path.join(cwd, 'last_simulation_name.dat'),'r').read().strip())

## TODO generate time-domain damped oscillators & test
## TODO decimate the freq axis, and test 

## Load time, real field and imaginary field
x, Eabs, Ephase = np.loadtxt(last_simulation_name+'_timedomain.dat', usecols=list(range(3)), unpack=True)
#y = Eabs * np.exp(1j*Ephase) * np.exp(-x/1e-9)


omega0, gamma, ampli = 3e9*2*np.pi, 1e9, 1e6
y = (np.sign(x)/2+.5) * np.sin(x*omega0)*ampli * np.exp(-x*gamma/2)  ## damped oscillator

#y = np.sin(2*np.pi*x*3e9) * np.exp(-x*3e8) * 1e6
## if oscillator's time constant is 1/100 times angular frequency, its quality is 100pi
def lorentz(omega, omega0, gamma, ampli):
    return 1/(2*np.pi) * ampli * omega0**2 / (omega0**2 - omega**2 + 1j*omega*gamma) 




## Convert to polar notation and save the spectrum
#meep_utils.loadtxt(fname=model.simulation_name+"_freqdomain.dat", X=zip(freq, np.abs(yf), meep_utils.get_phase(yf)), fmt="%.6e",
        #header=model.parameterstring + meep_utils.sim_param_string(sim_param) + "#x-column _frequency [Hz]\n#column ampli\n#column phase\n")

## Plot time-domain data
import matplotlib.pyplot as plt
plt.figure(figsize=(10,10))
#plt.plot(x, np.abs(y), color="k", label=u"$|y|$", ls='-')       
plt.plot(x, np.abs(y), color="k", label=u"y", ls='-')       
plt.xlabel(u"time [s]"); plt.ylabel(u"amplitude"); plt.yscale('log'); plt.grid()
plt.legend(prop={'size':10}, loc='upper right').draw_frame(False)
plt.savefig("td.png", bbox_inches='tight')

## 1D FFT with cropping for useful frequencies
plt.figure(figsize=(15,10))
freq    = np.fft.fftfreq(len(x), d=(x[1]-x[0]))         # calculate the frequency axis with proper spacing
yf      = np.fft.fft(y, axis=0) / len(x) * 2*np.pi      # calculate the FFT values
freq    = np.fft.fftshift(freq)                         # ensures the frequency axis is a growing function
yf      = np.fft.fftshift(yf) / np.exp(1j*2*np.pi*freq * x[0])   # dtto, and corrects the phase for the case when x[0] != 0
truncated = np.logical_and(freq>0, freq<maxftmp)         # (optional) get the frequency range
(yf, freq) = map(lambda array: array[truncated], (yf, freq))    # (optional) truncate the data points
plt.plot(freq, np.abs(yf), color="#FF8800", label=u"FFT", ls='-')                  # (optional) plot amplitude

## Convert to polar notation and save the spectrum
#meep_utils.savetxt(fname=model.simulation_name+"_freqdomain.dat", X=zip(freq, np.abs(yf), meep_utils.get_phase(yf)), fmt="%.6e",
        #header=model.parameterstring + meep_utils.sim_param_string(sim_param) + "#x-column _frequency [Hz]\n#column ampli\n#column phase\n")


## Harminv
## TODO switch to harminv_wrapper.py instead
#hi = meep_utils.harminv(x, y, amplitude_prescaling=1e6)
#oscillator_count = len(hi['frequency'])
import harminv_wrapper
hi = harminv_wrapper.harminv(x, y, amplitude_prescaling=1e0)

if type(hi['frequency']) == np.float64:
    hi['frequency'], hi['decay'], hi['amplitude'] = np.array([hi['frequency']]), np.array([hi['decay']]), np.array([hi['amplitude']])
oscillator_count = len(hi['frequency'])
if oscillator_count > 0:
    plt.scatter(np.abs(hi['frequency']), hi['amplitude'], c=hi['phase'], s=np.abs(hi['quality'])/50 + 2, cmap=plt.cm.hsv, alpha=.2)
    #for osc in zip(np.abs(hi['frequency']), hi['amplitude'], np.abs(hi['quality'])):
        #print 'f=%g A=%g q=%g' % osc

freq_fine = np.linspace(0, np.max(freq)*1, 10000)

plt.plot(freq_fine, 
        np.abs(lorentz(freq_fine*2*np.pi, omega0, gamma, ampli)), 
        color="#88FF00", label=u"Aosc", ls='-')      # (optional) plot amplitude

sumosc = np.zeros_like(freq_fine)*1j
print "Harminv frequencies", hi['frequency']
for osc in range(oscillator_count):
    osc_y = lorentz(freq_fine,   hi['frequency'][osc], hi['decay'][osc], hi['amplitude'][osc])
    plt.plot(freq_fine, np.abs(osc_y), color="#0088FF", label=u"", ls='-', alpha=.3)      # (optional) plot amplitude
    sumosc += osc_y 
plt.plot(freq_fine, np.abs(sumosc), color="#0088FF", label=u"$\\Sigma$ osc", ls='-')      # (optional) plot amplitude

analytic_modes = {}
from scipy.special import jnyn_zeros

#radius=33.774e-3
#height=122.36e-3
#freq_correction =  1       # (1. - 150e6/3.2e9)   ## optionally compensate for the few percent error in frequency (introduced by discretisation in FDTD) 
#for p in [0,1,2]:
    #for n in range(4):
        #S = " "*p
        #for m,B in enumerate(jnyn_zeros(n, 5)[0]):
            #analytic_modes[freq_correction * c/(2*np.pi) * np.sqrt((B/(radius))**2 + (p*np.pi/height)**2) ] = ("%s$TM_{%d%d%d}$%s" % (S,n,m+1,p,S))
        #for m,B in enumerate(jnyn_zeros(n, 5)[1]):
            #if p>0: ## TExx0 modes can not exist [Pozar]
                #analytic_modes[freq_correction * c/(2*np.pi) * np.sqrt((B/(radius))**2 + (p*np.pi/height)**2) ] = ("   %s$TE_{%d%d%d}$%s   " % (S,n,m+1,p,S))
#print analytic_modes
#meep_utils.annotate_frequency_axis(analytic_modes, label_position_y=1, arrow_length=10, log_y=True)

## Finish the plot + save 
plt.xlabel(u"frequency [Hz]")
plt.ylabel(u"amplitude excited by a pulse") 
plt.xlim((0, maxftmp))
#plt.ylim((1e1, 1e9))
plt.yscale('log')
plt.grid()
plt.legend(prop={'size':10}, loc='upper right').draw_frame(False)
plt.savefig("OUT.png", bbox_inches='tight')


#meep_utils.savetxt(freq=freq, s11=s11, s12=s12, model=model)
#with open("./last_simulation_name.dat", "w") as outfile: outfile.write(model.simulation_name) 
#import effparam        # process effective parameters for metamaterials

meep.all_wait()         # Wait until all file operations are finished








#def lorentz(x, f, d, q, A, p, err):
    #return A*abs(q)/((f/(2*np.pi))**2 + abs(q/2)*(x-f)**2) / (np.pi*2)          * 100
    #return A/(d + (x-f)**2) 
    #return A*abs(q)/(1 + abs(q*np.pi*2)*(x-abs(f))**2) / (np.pi*2) * 1e16 # XXX

