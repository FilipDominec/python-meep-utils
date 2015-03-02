#!/usr/bin/env python
#coding:utf8 
"""
Here you can find various functions and classes that facilitate the work with python-meep.
I believe some of these functions ought to be implemented in the meep module.
Filip Dominec 2012-2015

TODOs optional:
    * allow simple switching between 3-D and 2-D simulation
            1) requires to send the field along X direction, not Z
            2) choosing between 2-D field polarisations needs rotating source/monitor 
                i.e. (Ez,Hy) for TM, and  (Ey,Hz) for TE
            3) monitor plane must be *fast*, so it should switch to 1-D averaging at startup
    * try to make meep use true speed of light 
            1) multiply permittivity by epsilon0, and 
            2) introduce (flat) permeability mu0
            3) in forw/back wave separation, use true vacuum impedance
            4) remove all divisions/multiplications by speed of light (hooray!!)
            5) adjust Courant factor
    * verify that everything works also in the GNU Screen session (without xserver)
      ...   matplotlib.use('Agg') ## Enable plotting even in the GNU screen session
    * think over a very generic class of 
      SliceAnalysis(Slices=[], TimeDomainFunction=None, FreqDomainFunction=None, plot_timedomain=True, plot_freqdomain=True)
      instead of the MonitorPlane
"""
import os, os.path, sys, subprocess, time
import numpy as np
from scipy.constants import c, epsilon_0, mu_0

import matplotlib
matplotlib.use('Agg') # may help against the "tkinter.TclError: bad screen distance" error
import matplotlib.pyplot as plt

import meep_mpi as meep
#import meep

## === User interaction and convenience routines ===
def process_param(args):#{{{                
    """ Parse command-line parameters, and separate them into two groups.

    Some of them control the simulation (`sim_param'), but all remaining will be passed to 
    the model (`model_param')
    """
    sim_param = {   'frequency_domain':False,
                    'frequency':       None,
                    'MaxIter':         5000,
                    'MaxTol':          1e-2,
                    'BiCGStab':        8 }      ## BiCGStab order of 8 proved to have best performance
    model_param = {}
    for namevalue in args: ## first filter out those parameters that are specific for the simulation, rather than the model
        name, value = namevalue.split("=")
        if name == "frequency": 
            sim_param['frequency']          = phys_to_float(value)
            sim_param['frequency_domain']   = True
        elif name == "maxtol":      sim_param['MaxTol']     = phys_to_float(value)
        elif name == "maxiter":     sim_param['MaxIter']    = int(value)
        elif name == "bicgstab":    sim_param['BiCGStab']   = int(value)
        elif name == "Kx":          sim_param['Kx']         = phys_to_float(value)
        elif name == "Ky":          sim_param['Ky']         = phys_to_float(value)
        elif name == "Kz":          sim_param['Kz']         = phys_to_float(value)
        else:           ## all other parameters will be passed to the model:
            model_param[name] = phys_to_float(value)
            if type(model_param[name]) == str:
                meep.master_printf("  <str>   %s%s = %s\n" % (name, " "*max(10-len(name), 0), model_param[name]))
            else:
                meep.master_printf("  <float> %s%s = %.3e\n" % (name, " "*max(10-len(name), 0), model_param[name]))
    return sim_param, model_param
#}}}
def phys_to_float(s):#{{{
    """
    Float() that also recognizes the short SI prefixes. Returns string if value is not a number.
    >>> phys_to_float('0.0121')
    0.0121
    >>> phys_to_float('12.1m')
    0.0121
    >>> phys_to_float('121e-4')
    0.0121
    >>> phys_to_float('abcd')
    'abcd'
    """
    prefixes = {'z':1e-21, 'a':1e-18, 'f':1e-15, 'p':1e-12, 'n':1e-9, 'u':1e-6, 'm':1e-3, 
            'c':1e-2, 'd':1e-1, 'k':1e3, 'M':1e6, 'G':1e9, 'T':1e12, 'P':1e15, 'E':1e18, 'Y':1e21}
    try:
        if s[-1] in '.0123456789':       
            return float(s)
        elif s[-1] in prefixes.keys():                                  
            return float(s[:-1]) * prefixes[s[-1]]
        else:
            return s
    except ValueError:
        return s#}}}
class Timer():#{{{ 
    """
    Prints total estimated time of computation, and the simulation progress 
    """
    def __init__(self, simtime):
        self.starttime = time.time()
        self.simtime = simtime
        meep.master_printf("\tSimulation time: %e [s] = %e time units\n" % (simtime, simtime*c))
        self.reporttimes = [.001, .01, 0.03] + [t/10. for t in range(1,10)] + [2.]
    def get_time(self):
        return time.time()-self.starttime
    def print_progress(self, now): 
        if now/self.simtime > self.reporttimes[0]:
            meep.master_printf("Progress %.2f of expected total %d s\n" % \
                    (now / self.simtime, (self.simtime / now * self.get_time())))
            self.reporttimes[0:1] = []
#}}}
def notify(title, run_time=None):#{{{
    """
    Shows a bubble with notification that your results are about to be ready!
    Requires python-notify installed, otherwise just quits

    Note: you may also get similar results with libnotify-bin from a bash script:
        notify-send -t 3000 -i "face-glasses" "something happened"
    """
    if meep.my_rank() != 0: return
    try: 
        if run_time: timestring = "in %d s" % int(run_time)
        else: timestring = ""
        import pynotify
        pynotify.init("image")
        n = pynotify.Notification("MEEP simulation finished %s" % (timestring), title, "face-glasses")
        n.show()
    except:
        pass
#}}}
def last_simulation_name(argindex=1): #{{{ 
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
def find_maxima(x, y, minimum_value=0):#{{{
    """ 
    Returns the x points where 
    1) y has a local maximum (i. e. dx/dy goes negative) AND 
    2) where y is above minimum_value 
    """
    d = y[1:-1] - y[0:-2]   ## naïve first derivative
    maxima = x[1:][np.sign(d[0:-2])-np.sign(d[1:-1]) + np.sign(y[2:-2]-minimum_value)==3]
    return maxima 
#}}}


## === Structure definition and setup ===
## Define the simulated models as a class (much simpler than handling callbacks manually)
class AbstractMeepModel(meep.Callback): 
    def __init__(self):#{{{
        meep.Callback.__init__(self)
        self.double_vec = None          # (callback function to be redirected to the desired function)
        self.return_value = True  
        #}}}
    def eps(self, r):#{{{
        """ Scans through materials and returns the high-frequency part of permittivity for the first in the list. 
        Be careful when materials overlap - their polarizabilities still sum up, making a nonrealistic result. 
        """
        ## TODO #CSGspeedup rewrite to:  sum([mat.eps for mat in self.materials if mat.where(r)]) or 1 
        ## TODO rewrite the geometric functions to lambdas and test the speed!
        for mat in self.materials:
            if mat.where(r): return mat.eps         ## TODO: #realCproject will MEEP use the natural speed of light 3e8 m/s, if this is multiplied by eps0?
        else: return 1.                             ## TODO: and here too   (... this also needs that similar approach is used to define mu=mu0 everywhere)
        #}}}
    def register_local(self, param, val):#{{{
        """ Adds a parameter as an attribute of the model (either number or float). 
        If the parameter value differs its default one, adds it also to the simulation name"""
        setattr(self, param, val)

        nondefault = self.named_param_defaults.get(param, None) != val ## XXX
        if param not in self.named_param_defaults.keys(): infostring = "(set in code)"
        elif nondefault: infostring = "" 
        else: infostring = "(default)" 

        ## prepare the parameter to be added into name (if not conversible to float, add it as a string)
        try: 
            if nondefault: self.simulation_name += ("_%s=%.3e") % (param, float(val))
            self.parameterstring += "#param %s,%.4e\n" % (param, val)
            #meep.master_printf("  <float> %s%s = %.3e %s\n" % (param, " "*max(10-len(param), 0), val, infostring))
        except ValueError: 
            if nondefault: self.simulation_name += ("_%s=%s") % (param, val)
            self.parameterstring += "#param %s,%s\n" % (param, val)
            #meep.master_printf("  <str>   %s%s = %s %s\n" % (param, " "*max(10-len(param), 0), val, infostring))
        #}}}
    def register_locals(self, params):#{{{
        """ Scans through the parameters and calls register_local() for each """
        import inspect
        a = inspect.getargspec(self.__init__)
        self.named_param_defaults = dict(zip(a.args[-len(a.defaults):], a.defaults))

        self.parameterstring = ""
        ## First look up for the "VIP" parameters that should come first in the name:
        preferred_params = ['resolution', 'comment', 'frequency', 'simtime']
        for param in preferred_params:
            if params.get(param): 
                val = params.get(param)
                self.register_local(param, val)
        ## Then add all remaining parameters of the model
        for (param, val) in params.iteritems():
            if param != 'self' and param not in preferred_params:
                self.register_local(param, val)
        #}}}
    def f_c(self):#{{{
        """ critical_frequency for FDTD stability """
        return c / np.pi/self.resolution/meep.use_Courant() 
        #}}}
    def build_polarizabilities(self, structure):#{{{
        """ 
        This is a helper to define the susceptibilities using the callback.
        It goes through all polarizabilities for all materials. 
        Applicable for time-domain simulation only, because dispersive model is not implemented for 
        frequency-domain simulation yet.
        """
        avail_cbs = [meep.DBL5, meep.DBL4, meep.DBL3, meep.DBL2, meep.DBL1,]
        avail_cb_setters = [meep.set_DBL5_Callback, meep.set_DBL4_Callback, meep.set_DBL3_Callback, 
                meep.set_DBL2_Callback, meep.set_DBL1_Callback,]
        for material in self.materials:
            aeps = analytic_eps(material, self.src_freq)
            meep.master_printf("Info\tAdding material: %s with %d oscillator(s); (eps @ %.2e Hz = %.1f+%.3fj)\n" % 
                    (material.name, len(material.pol), self.src_freq, aeps.real, aeps.imag))
            self.double_vec = material.where  ## redirect the double_vec() function callback
            for polariz in material.pol:
                if avail_cbs == []: 
                    meep.master_printf("Error: too many oscillators defined in total (>5)."+"Perhaps the simulation can be made with less?")
                    exit()
                next_cb, next_cb_setter = avail_cbs.pop(), avail_cb_setters.pop()
                self.return_value = polariz['sigma']
                next_cb_setter(self.__disown__())    
                if "lorentzian_susceptibility" in dir(meep):
                    ## for meep 1.2 or newer
                    structure.add_susceptibility(next_cb, meep.E_stuff, 
                            meep.lorentzian_susceptibility(polariz['omega']/c, polariz['gamma']/c))     ## todo #realCproject":
                    #else:todo: fix in python-meep
                        #print dir(meep)
                        #structure.add_susceptibility(next_cb, meep.E_stuff, 
                                #meep.drude_susceptibility(polariz['omega']/c, polariz['gamma']/c)) 
                else:
                    ## for meep 1.1.1 or older
                    structure.add_polarizability(next_cb, polariz['omega']/c, polariz['gamma']/c)       ## todo #realCproject":

            # TODO  fix crashing 
            #if 'chi2' in dir(material):
                #self.return_value = material.chi2
                #meep.set_CHI2_Callback(self.__disown__())    
                #structure.set_chi2(meep.CHI2)
            #if 'chi3' in dir(material):
                #self.return_value = material.chi3
                #meep.set_CHI3_Callback(self.__disown__())    
                #structure.set_chi3(meep.CHI3)

            # #structure.set_chi3(meep.Ex, None)
            # # Todo: can a vector or 4-D tensor be used?
            # # structure.set_chi3(meep.X, meep.CHI3) ... is wrong, because:
            # #  --->  NotImplementedError: Wrong number of arguments for overloaded function 'structure_orig_set_chi3'.
            #           #Possible C/C++ prototypes are:
            #             #set_chi3(meep::structure *,meep::component,meep::material_function &)
            #             #set_chi3(meep::structure *,meep::material_function &)
            #             #set_chi3(meep::structure *,double (*)(meep::vec const &)) XXX this is perhaps used XXX
        #}}}
    def fix_material_stability(self, material, f_c="Auto", verbose="false"):#{{{
        """ Little heuristics to make the time-domain simulation (mostly) stable

        First, all oscillators with too high frequencies are removed, and the nondispersive
        part of permittivity is increased. Second, the possible Drude term is detected and fixed
        so that metals work in low-resolution simulations as well.

        Note that this function is not guarranteed to give optimal results. Sometimes an oscillator 
        is so strong that it pulls permittivity negative above the critical frequency f_c, and 
        simulation goes unstable. Very often the models need to be adjusted for the simulation to
        be faster and more accurate.

        More information is on the authors website.
        """
        if f_c == "Auto": f_c = self.f_c()
        f_c_safe = f_c * 0.5

        ## Check and fix the first stability criterion
        for n, osc in enumerate(material.pol[:]):
            if osc['omega'] > f_c_safe: 
                material.eps += osc['sigma']
                material.pol.remove(osc)
                if verbose: 
                    meep.master_printf("Removing oscillator #%d at too a high frequency (%.2e) from material: %s\n" % \
                            (n+1, osc['omega'], material.name))

        ## Find possible Drude terms and fix them if they break the 2nd stability criterion
        for osc in material.pol:
            if osc['omega'] < 1:      ## detects Drude term, TODO make more general (or even define new parameters for Drude term)
                ## Retrieve the properties of the Drude term
                gamma    = osc['gamma'] * 2*np.pi  # angular scattering frequency 

                ## Real part of permittivity is roughly constant under the scattering frequency, `gamma', 
                ## and grows above it. To push permittivity above zero at f_c, gamma must be lower than f_c
                ## but not too much. One half is reasonable.
                max_gamma = f_c_safe
                if osc['gamma'] > max_gamma: 
                    osc['gamma'] = max_gamma

                ## For a typical Drude conductive material, permittivity is such a big negative number at low frequencies,
                ## so pushing it above 1 at high frequencies does not make any appreciable error.
                ## To be on the safe side, we ensure eps' > 1 even at f_c*0.5.
                eps_at_fc = np.real(analytic_eps(material, freq=f_c * 0.5))
                if eps_at_fc < 1:
                    material.eps += 1 - eps_at_fc
                    if verbose: 
                        meep.master_printf(("Increasing high-frequency permittivity by %.1f to "+ \
                                "make stable the material: %s\n") % (1 - eps_at_fc, material.name))
        #}}}
    def test_materials(self, verbose="false"):#{{{
        """ 
        1) Verify the material definition will not introduce instabilities in FDTD
        2) Call the where() function for each material, in order to make sure there are no errors
        (SWIG callback does not report where the error occured, it just crashes) 
        """
        f_c = self.f_c()
        for n, material in enumerate(self.materials): 
            ## Check the first stability criterion (that MEEP checks, but in a wrong way)
            for osc in material.pol:
                if osc['omega'] > f_c: 
                    meep.master_printf("\n\tWARNING: the oscillator %d in material `%s' defined above the critical frequency (%g)." % (n, material.name, f_c))
                    meep.master_printf("\n\t         It may help to run fix_material_stability for this material, or edit it manually.\n\n")

            ## Check the 2nd stability criterion for all remaining oscillators
            eps_fc      = analytic_eps(self.materials[0], f_c)
            eps_minimum = meep.use_Courant()**2 * 3 # for three dimensions (2-D simulations are safer as they multiply by 2 only, but it is ignored here)
            if (eps_fc.real < eps_minimum):
                meep.master_printf("\n\tWARNING: at the critical frequency %f, real permittivity of %s is below the criterion for stability (%f < %f)"
                        % (f_c, material.name, eps_fc, eps_minimum))

            ## Test whether the `where()' function 
            if callable(material.where):    
                for x in np.linspace(-self.size_x/2, self.size_x/2, 15):
                    for y in np.linspace(-self.size_y/2, self.size_y/2, 15):
                        if self.size_z:     # 3D case
                            for z in np.linspace(-self.size_z/2, self.size_z/2, 15): material.where(meep.vec(x, y, z))
                        else:               # 2D case
                            material.where(meep.vec(x, y))
            else:       
                meep.master_printf("\n\tWARNING: `where' parameter is not a function, material not used: %s\n\n" % (material.name))
                self.materials.remove(material)
        #}}}

## Geometrical primitives to help defining the geometry
def in_xslab(r,cx,d):#{{{
    return (abs(r.x()-cx) < d/2)
def in_yslab(r,cy,d):
    return (abs(r.y()-cy) < d/2)
def in_zslab(r,cz,d):
    return (abs(r.z()-cz) < d/2)
def in_xcyl(r,cy,cz,rad):
    return ((r.y()-cy)**2+(r.z()-cz)**2) < rad**2
def in_ycyl(r,cx,cz,rad):
    return ((r.x()-cx)**2+(r.z()-cz)**2) < rad**2
def in_zcyl(r,cx,cy,rad):
    return ((r.x()-cx)**2+(r.y()-cy)**2) < rad**2
def in_sphere(r,cx,cy,cz,rad):
    return ((cx-r.x())**2 + (cy-r.y())**2 + (cz-r.z())**2)**.5 < rad
def in_ellipsoid(r,cx,cy,cz,rad,ex):
    xd, yd, zd = (cx-r.x()), (cy-r.y()), (cz-r.z())
    return ((xd+yd)**2/2.*ex**2 + (xd-yd)**2/2./ex**2 + zd**2)**.5 < rad
#}}}

## Note: An example of how to define a custom source - put this code into your simulation script 
#{{{
"""
## Replace the f.add_volume_source(meep.Ex, srctype, srcvolume) line with following:
## Option for a custom source (e.g. exciting some waveguide mode)
class SrcAmplitudeFactor(meep.Callback): 
    ## The source amplitude is complex -> phase factor modifies its direction
    ## todo: implement in MEEP: we should define an AmplitudeVolume() object and reuse it for monitors later
    def __init__(self, Kx=0, Ky=0): 
        meep.Callback.__init__(self)
        (self.Kx, self.Ky) = Kx, Ky
    def complex_vec(self, vec):   ## Note: the 'vec' coordinates are _relative_ to the source center
        # (oblique) plane wave source:
        #return np.exp(-1j*(self.Kx*vec.x() + self.Ky*vec.y()))
        # (oblique) Gaussian beam source:
        return np.exp(-1j*(self.Kx*vec.x() + self.Ky*vec.y()) - (vec.x()/100e-6)**2 - (vec.y()/100e-6)**2) 
af = SrcAmplitudeFactor(Kx=model.Kx, Ky=model.Ky)
meep.set_AMPL_Callback(af.__disown__())
f.add_volume_source(meep.Ex, srctype, srcvolume, meep.AMPL)

## BTW could one not use an anonymous object here? How does it cope with this? meep.set_AMPL_Callback(af.__disown__())
## o = type('AmplitudeFactor', (meep.Callback,), { "complex_vec": lambda xxx})
## ...
"""
#}}}


## === Initialization of the materials, structure and whole simulation ===
def permittivity2conductivity(complex_eps, freq):#{{{
    """
    Enables to use the same dispersive materials for time- and frequency-domain simulation

    Complex permittivity can express also the conductivity of the sample (in the same
    manner as dielectric losses) with the corresponding relation:
        complex_eps = real_eps - 1j conductivity / (frequency * 2*pi * epsilon_0)
    Therefore it should be inverted for classic D-conductivity:
        conductivity = -
    In order to simulate any lossy medium with the freq-domain solver, we invert this relation
    to obtain a (nondispersive) conductivity for one frequency. But it does not give the same results
    as time-domain simulation.

        What we know: 
           (we know that     c**.5 * np.pi  = 54395 ) 
        function of c, f, 2pi, eps0, eps.im/eps.r
        should give dimension 1  to feed unitless meep
        should be proportional to eps.im/eps.r
        should give ca. 50000 for omega = 2pi * 800 GHz and eps.im/eps.r=0.02
           => should give 2.5e6 for (eps.im/eps.r=1)
        should be proportional to frequency omega
            => should give 5e-7 for omega = 1 and  (eps.im/eps.r=1)
        already was pre-divided for meep by c = 3e8  (which counts here)
            => in real life it shall be  3e-6 * 3e8 = 148
        should be proportional to epsilon0 [C/Vsm], which is similar in units to conductivity
            => if epsilon0 was 1, it would be 1.7e13 -> roughly c**2
    """
    ## TODO resolve the Bulgarian constant when running freq-domain simulation
    # return complex_eps.imag * freq * 2*np.pi * epsilon_0 * complex_eps.real  ## orig. idea
    # return complex_eps.imag * freq * 2*np.pi * epsilon_0 * complex_eps.real ## also wrong
    # return complex_eps.imag / complex_eps.real * 2*np.pi * c
    #return complex_eps.imag / complex_eps.real * 6.28*freq * 8.85e-12 * c
    magic_constant = 1.65e13       ## A. K. A. bulgarian constant...
    return complex_eps.imag / complex_eps.real * 2 * np.pi * freq * epsilon_0 * c**.5 * np.pi
#}}}
def analytic_eps(mat, freq):#{{{
    complex_eps = np.ones_like(freq)*(mat.eps+0j)
    #omega = freq * 2*np.pi
    for polariz in mat.pol:
        complex_eps += polariz['sigma'] * polariz['omega']**2 / (polariz['omega']**2 - freq**2 - 1j*freq*polariz['gamma']) 
        #complex_eps += polariz['sigma'] * polariz['omega']**2 / (polariz['omega']**2 - omega**2 - 1j*omega*polariz['gamma']) 
    return complex_eps 
#}}}
class MyHiFreqPermittivity(meep.Callback):#{{{      %% TODO rename to Permittivity_callback
    def __init__(self, model, frequency):
        meep.Callback.__init__(self)
        self.model = model
        self.frequency = frequency
    def double_vec(self, r):
        for material in self.model.materials:
            if material.where(r):
                return analytic_eps(material, self.frequency).real
        else: return 1
#}}}
class MyConductivity(meep.Callback):#{{{            %% TODO rename to Conductivity_callback
    def __init__(self, model, frequency):
        meep.Callback.__init__(self)
        self.model = model
        self.frequency = frequency
    def double_vec(self, r):
        for material in self.model.materials:
            if material.where(r):
                return permittivity2conductivity(analytic_eps(material, self.frequency), self.frequency)
        else: return 0
#}}}

def annotate_frequency_axis(mark_freq, label_position_y=1, arrow_length=3, log_y=False):#{{{
    """
    """
    import matplotlib.pyplot as plt
    if type(mark_freq) in (float,int): mark_freq=list(mark_freq,)
    if type(mark_freq) in (list,tuple): mark_freq=dict(zip(mark_freq,['' for mf in mark_freq]))
    for mfreq, mfreqtxt in mark_freq.items(): 
        label_y2 = label_position_y
        while (mfreqtxt[0:1]==' ' and mfreqtxt[-2:-1]==' '): 
            label_y2 = label_y2*2 if log_y else label_y2+1
            mfreqtxt=mfreqtxt[1:-1]; 
        bboxprops   = dict(boxstyle='round, pad=.15', fc='white', alpha=1, lw=0)
        arrowprops  = dict(arrowstyle=('->', '-|>', 'simple', 'fancy')[0], connectionstyle = 'arc3,rad=0', lw=1, ec='k', fc='w')
        plt.annotate(mfreqtxt,                    
                xy      = (mfreq, label_y2),    xycoords  ='data',
                xytext  = (mfreq, label_y2*arrow_length if log_y else label_y+arrow_length),  textcoords='data',        # (delete this if text without arrow is used)
                ha='center', va='bottom', size=15, color='k',
                bbox        = bboxprops,        # comment out to disable bounding box
                arrowprops  = arrowprops,       # comment out to disable arrow
                )
#}}}
def plot_eps(to_plot, filename="epsilon.png", plot_conductivity=True, freq_range=(1e10, 1e18), mark_freq=[], draw_instability_area=None):#{{{
    """ Plots complex permittivity of the materials to a PNG file

    Accepts list of materials
    """

    matplotlib.rc('text', usetex=True)
    matplotlib.rc('text.latex', preamble = '') # \usepackage{lmodern}, 

    #for material in list(to_plot): ## autoscale x axis?
        #for pol in material.pol:
            #if freq_range[1] < pol['omega']: freq_range[1] = pol['omega']*2

    frequency = 10**np.arange(np.log10(freq_range[0]), np.log10(freq_range[1]), .01)

    plt.figure(figsize=(7,6))
    #colors = ['#000000', '#004400', '#003366', '#000088', '#440077', '#661100', 
              #'#aa8800', '#00aa00', '#0099dd', '#0000EE', '#2200DD', '#aa0000']
    colors = ['#000000', '#004400', '#003366', '#000088', '#440077', '#661100', 
              '#aa8800', '#0044dd', '#00bb00', '#aaaa00', '#bb6600', '#dd0000']

    subplotnumber = 2 if plot_conductivity else 1


    for material in list(to_plot):
        if colors: color = colors.pop()
        else: color = 'black'
        label = getattr(material, 'shortname', material.name)
        plt.subplot(subplotnumber,1,1)

        eps = np.conj(analytic_eps(material, frequency)) ## FIXME eps should be computed as conjugated by default

        plt.subplot(subplotnumber,1,1)
        plt.plot(frequency, np.real(eps), color=color, label=material.name, ls='-') #  
        plt.plot(frequency, np.imag(eps), color=color, label='', ls='--') # 
        #R = abs((1-eps**.5)/(1+eps**.5))**2     ## Intensity reflectivity

        if plot_conductivity:
            plt.subplot(subplotnumber,1,2)
            label = ""
            omega = 2*np.pi*frequency
            cond = eps * omega * epsilon_0 * 1j ## for positive frequency convention
            #cond = eps * omega * epsilon_0 / 1j ## for negative frequency convention

            plt.plot(frequency, np.real(cond), color=color, label=material.name, ls='-')
            plt.plot(frequency, np.imag(cond), color=color, label="#", ls='--')


            #plt.plot(frequency, np.real(1./cond)*1e12, color=color, lw=1.5, label=("$\\rho^{'}\cdot 10^{12}$"), ls=':') ## (real) resistivity

            # XXX temporary
            #gamma = 1.5e+13
            #f_p = 2.07153080589e+14
            #omega_p =  1.3015811923e+15
            #omega_p = 1.3015811923e+15 # / (2*np.pi)**.5

            ## Low-frequency limits for pseudo-Drude model
            #plt.plot(frequency, np.ones_like(omega) * omega_p**2 * epsilon_0 / gamma, color='k', label=label, ls='-', lw=.3)
            #plt.plot(frequency, omega * (-epsilon_0 + omega_p**2 * epsilon_0 / gamma**2), color='k', label=label, ls='--', lw=.3)
            ## High-frequency limits for pseudo-Drude model
            #plt.plot(frequency, omega**-2 * (omega_p**2 * epsilon_0 * gamma), color='g', label=label, ls='-', lw=.3)
            #plt.plot(frequency, omega * -epsilon_0            , color='g', label=label, ls='--', lw=.3)
                                                    #^??????^/(2*np.pi)
            # XXX                   #  ^^ ?????????? ^

            ## TODO check this and perhaps remove:
            #plt.subplot(subplotnumber,1,3)
            #plt.plot(frequency, -np.real(cond), color=color, label=label, ls='-')
            #plt.plot(frequency, -np.imag(cond), color=color, label="", ls='--')
            #
            #
            #plt.ylabel(u"negative valued")
            #plt.yscale('log'); 
            #plt.xscale('log'); plt.legend(); plt.grid(True)



    ## Annotate frequencies and finish the graph 
    plt.subplot(subplotnumber,1,1)
    plt.legend(prop={'size':10}) 
    plt.xlabel(u"frequency $f$ [Hz]") 
    plt.ylabel(u"relative permittivity $\\varepsilon_r$")
    plt.xscale('log'); plt.grid(True)
    ylim = (-1e7, 1e6); plt.ylim(ylim); plt.yscale('symlog')
    annotate_frequency_axis(mark_freq, log_y=True, arrow_length=50) # TODO , print_freq=True
    if draw_instability_area:
        plt.gca().add_patch(plt.Rectangle((draw_instability_area[0], ylim[0]), 1e20, draw_instability_area[1]-ylim[0], color='#dddddd'))

    if plot_conductivity:
        plt.subplot(subplotnumber,1,2)
        plt.ylabel(u"conductivity $\\sigma$")
        plt.xscale('log'); plt.grid(True)
        plt.yscale('symlog'); plt.ylim((-1e12, 1e12)); 
        plt.xlabel(u"frequency $f$ [Hz]") 

    plt.xlabel(u"frequency $f$ [Hz]") 
    plt.savefig(filename, bbox_inches='tight')
    plt.savefig(filename+".pdf", bbox_inches='tight')
#}}}

def init_structure(model, volume, sim_param, pml_axes):#{{{
    """
    This routine wraps the usual tasks needed to set up a realistic simulation with meep.

    model, volume and sim_param are objects that need to be passed from the main simulation

    pml_axes may be selected from those: None, meep.X, meep.XY, meep.Y or meep.Z, "All" 
    """
    def init_perfectly_matched_layers():
        if pml_axes == "All" or pml_axes == "all":
            perfectly_matched_layers=meep.pml(model.pml_thickness)
            s = meep.structure(volume, meep.EPS, perfectly_matched_layers, meep.identity())
        elif pml_axes == None or pml_axes == "none" or pml_axes == "None":
            s = meep.structure(volume, meep.EPS)
        else:
            perfectly_matched_layers=meep.pml(model.pml_thickness, pml_axes)
            s = meep.structure(volume, meep.EPS, perfectly_matched_layers, meep.identity())
        return s

    if not sim_param['frequency_domain']:
        meep.master_printf("== Time domain structure setup ==\n")
        ## Define each polarizability by redirecting the callback to the corresponding "where_material" function
        ## Define the frequency-independent epsilon for all materials (needed here, before defining s, or unstable)
        model.double_vec = model.eps
        meep.set_EPS_Callback(model.__disown__())
        structure = init_perfectly_matched_layers()

        ## Add all the materials
        model.build_polarizabilities(structure)

        ## TODO 1. Add chi2, chi3 to AbstractMeepModel (done?), 
        ## TODO 2. move  the callback-like definition from build_polarizabilities() here - similarly as with EPS
    else:
        meep.master_printf("== Frequency domain structure setup (for frequency of %g Hz) ==\n" % sim_param['frequency'])
        if (model.Kx!=0 or model.Ky!=0): print "Warning: frequency-domain solver may be broken for nonperpendicular incidence"
        ## Frequency-domain simulation in meep does not accept dispersive materials yet. We must redefine each material by 
        ## using the nondispersive permittivity and the nondispersive conductivity 
        ## (both calculated from polarizabilities at given frequency)

        ## Define the frequency-independent epsilon for all materials (has to be done _before_ defining s, otherwise unstable)
        my_eps = MyHiFreqPermittivity(model, sim_param['frequency'])
        meep.set_EPS_Callback(my_eps.__disown__())
        structure = init_perfectly_matched_layers()
        # TODO display a warning about freq-domain neglecting chi2, chi3 (if defined and nonzero in some of materials)

        ## Create callback to set nondispersive conductivity (depends on material polarizabilities and frequency)
        mc = MyConductivity(model, sim_param['frequency'])
        meep.set_COND_Callback(mc.__disown__())
        structure.set_conductivity(meep.E_stuff, meep.COND)  ## only "E_stuff" worked here for me

    return structure
#}}}


## === Results post-processing and export ===
## Call filter-diagonalisation method to analyze time-domain data  
## (TODO not using MEEP functions --> shall be synced against harminv_wrapper.py and removed from here)
def run_bash(cmd, anyprocess=False): #{{{
    if meep.my_rank() == 0 or anyprocess:
        meep.master_printf("CMD: "  + cmd+ "\n")  # (diagnostics)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        out = p.stdout.read().strip()
        return out
#}}}
def harminv(x, y, d=100, f=30, amplitude_prescaling=1):#{{{
    """
    This algorithm may be very sensitive to input amplitude etc.

    Suggested visualisation:
    hi = meep_utils.harminv(x,y)
    plt.scatter(hi['frequency'], hi['amplitude'], c=hi['phase'], s=np.abs(hi['quality'])/20 + 2, cmap=plt.cm.hsv, alpha=.3)
    """

    # TODO: try-except to catch when harminv returns no data (this happens)
    with open('hitest.dat', 'w') as outfile: 
        outfile.write("#t[s]\t E(t)\n")
        np.savetxt(outfile, zip(x/2, y * amplitude_prescaling), fmt="%.8e")
    import subprocess
    dt = x[1]-x[0]
    subprocess.Popen('harminv 0-%f -t %g -d %g -f %g < hitest.dat > hiout.dat' % (1/dt/10, dt, d, f), shell=True, stdout=subprocess.PIPE).stdout.read().strip()
    try:
        (mf, md, mQ, mA, mp, merr) = np.loadtxt('hiout.dat', usecols=list(range(6)), unpack=True, delimiter=', ', skiprows=1)
    except:
        meep.master_printf("\nWARNING: Harminv detected no resonances.\n\n")
        (mf, md, mQ, mA, mp, merr) = [np.array([]) for _ in range(6)]

    return {'frequency':np.abs(mf*2), 'decay':md, 'quality':mQ, 'amplitude':mA/amplitude_prescaling, 'phase':mp, 'error':merr}
#}}}

## Saving and loading data (not dependent on MEEP functions, but better if ran by the 1st process only)
def sim_param_string(sim_param):#{{{
    if sim_param['frequency_domain']:
        output = "#param frequency_domain,True\n#param SolverTol,%d\n#param SolverMaxIter,%d\n#param SolverBiCGStab,%d\n" % \
                (sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab'])
    else:
        output = "#param frequency_domain,False\n"
    if sim_param.get('Kx') != None: output += "#param Kx,%.3f\n" % sim_param.get('Kx')
    if sim_param.get('Ky') != None: output += "#param Ky,%.3f\n" % sim_param.get('Ky')
    if sim_param.get('Kz') != None: output += "#param Kz,%.3f\n" % sim_param.get('Kz')
    return output
#}}}
def savetxt(fname, X, header, **kwargs):#{{{ 
    """
    Its use is for older versions of the library that do not accept the `header' parameter
    """
    #FIXME - This function will be superseded by numpy.savetxt
    with open(fname, "w") as outfile: 
        outfile.write(header)
        np.savetxt(outfile, X, **kwargs)
#}}}
def loadtxt_params(filename):#{{{
    parameters  = {}
    with open(filename) as datafile:
        for line in datafile:
            if (line[0:1] in '0123456789') or ('column' in line.lower()): break    # end of parameter list
            key, value = line.replace(',', ' ').split()[-2:]
            try: value = float(value) ## Try to convert to float, if possible
            except: pass                ## otherwise keep as string
            parameters[key] = value
    return parameters
#}}}
def loadtxt_columns(filename):#{{{
    columns     = []
    with open(filename) as datafile:
        for line in datafile:
            if ('column' in line.lower()): columns.append(line.strip().split(' ', 1)[-1])
    return columns
#}}}

## Useful for making 3D snapshots of the fields
class Slice(): #{{{
    """ Exports a slice through the fourdimensional field information in space-time.
    If you specify time as a number, using the parameter `at_t', it takes a full 3D snapshot of the 
    field at that time. If the number exceeds the length of the simulation, the snapshot is taken at the end anyway.
    If you specify some the X, Y or Z coordinate using e.g. `at_x', it exports a time evolution in 
    the specified plane. More than one coordinate can be given, in which case the exported field has
    less dimensions, accordingly.

    The `at_x', `at_y', `at_z', `at_t' can also be of the <list> type with two elements, in which case 
    they delimit an interval in the volume or time and do not affect the dimension of the output.

    Some simulations can provide excessively long series at the time axis. You may set `min_timestep'
    to avoid exporting a data file in each time step, but in the defined.
    TODO: make it automatic, limiting the GIF by default to 5 seconds ~ 125 frames

    Output format can be changed by setting True one or more of these parameters:
        outputpng: directory of PNG images  (can produce too many files if `min_timestep' not properly set!)
            ## TODO: if more field components specified, export each of them separately, do not forget to add `.ex' to the name
        outputgif: one animated GIF
            ## TODO: if more field components specified, --dtto--
        outputhdf: one HDF file containing complex fields for all components
        outputvtk: real part of each component in separate VTK file (for mayavi2/paraview)

    -- Usage --
    One has to provide the model object with these attributes: size_x, size_y, size_z, simtime, simulation_name
    The components parameter may be one or more meep components, such as meep.Ex, meep.Ey or meep.Dielectric
    Also the fields object is compulsory.

    The name of the exported files is chosen automatically to avoid conflicts.

    -- Examples --
    >>> # Export dielectric structure:
    >>> meep_utils.Slice(model=model, field=f, components=(meep.Dielectric), at_t=0, outputhdf=True, name='EPS')]
    >>> # Export full-wave evolution in one X-Y plane
    >>> meep_utils.Slice(model=model, field=f, components=(meep.Ex, meep.Ey, meep.Ey), at_z=0, min_timestep=.1e-12, outputhdf=True, name="ZPLANE")]
    >>> # Export a final snapshot of the field, but restrict to one half of the simulation only
    >>> meep_utils.Slice(model=model, field=f, components=meep.Ex, at_t=model.simtime-.1e-12, at_x=[-np.inf, 10e-6], outputhdf=True, name="SNAP")
    """

    def isrange(self,tup): 
        """ Decides whether a space/time coordinate is a range (from-to), or just a single point """
        return tup[0] < tup[1]
    def __init__(self, model, field, components, 
            at_x=None, at_y=None, at_z=None, at_t=None, timebounds=None, min_timestep=0, 
            volume=None, pad=0, use_imag = False,
            outputdir=None, name='', outputpng=False, outputgif=False, outputhdf=False, outputvtk=True):
        def fix_xyzt_ranges(at, minlimit, maxlimit):#{{{
            """Returns a tuple: (from, to) """
            if type(at) in (tuple, list):       ## defined a range -> clip it to limits
                if at[0] > at[1]: at = [at[1], at[0]]
                if at[0] < minlimit: at[0] = minlimit    
                if at[1] > maxlimit: at[1] = maxlimit
                return at
            elif at == None:                    ## (default) no range or position defined
                return [minlimit, maxlimit]

            else:                               ## defined one position in space or in time
                return [float(at), float(at)]
        #}}}
        def generate_name(tuples, axis_names, component_name=None):#{{{
            nname=""
            ranges_count = 0
            for (tup, axis_name) in zip(tuples, axis_names):
                if tup and tup[0] == tup[1]: 
                    nname+="_%s%.3e" % (axis_name, tup[0])
            if component_name: name = name+"_"+component_name
            return nname
        #}}}

        self.outputpng, self.outputgif, self.outputhdf, self.outputvtk  = outputpng, outputgif, outputhdf, outputvtk
        self.field = field 
        self.components = (components,) if (type(components) not in (tuple, list)) else components
        self.real_or_imag = '.i' if use_imag else '.r'


        # (following assumes cartesian coords with volume.center_origin() called)
        at_x = fix_xyzt_ranges(at_x, -model.size_x/2, model.size_x/2)      
        at_y = fix_xyzt_ranges(at_y, -model.size_y/2, model.size_y/2)
        if model.size_z:
            at_z = fix_xyzt_ranges(at_z, -model.size_z/2, model.size_z/2)
        at_t = fix_xyzt_ranges(at_t, 0              , model.simtime )
        self.at_x, self.at_y, self.at_z, self.at_t = at_x, at_y, at_z, at_t
        self.min_timestep = min_timestep

        if model.size_z:
            self.volume = meep.volume(meep.vec(at_x[0], at_y[0], at_z[0]), meep.vec(at_x[1], at_y[1], at_z[1])) ## 3D
        else:
            self.volume = meep.volume(meep.vec(at_x[0], at_y[0]), meep.vec(at_x[1], at_y[1])) ## 2D
        
        #print at_x, at_y, at_z, at_t
        #print "Export dimension:", count_ranges((at_x, at_y, at_z, at_t))
        self.name = "%s_at%s" % (name if name else "Slice", 
                generate_name((self.at_x, self.at_y, self.at_z, self.at_t), 'xyzt'))

        if outputdir == None: outputdir = model.simulation_name
        if outputdir == None: outputdir = './'
        self.outputdir = outputdir

        self.images_number = 0
        self.last_slice_time = 0.

        meep.all_wait()
        if (meep.my_rank() == 0 and not os.path.exists(outputdir)): os.mkdir(outputdir)
        meep.all_wait()
        output_file_name = "%s.h5" % os.path.join(self.outputdir, self.name)
        if (meep.my_rank() == 0 and os.path.exists(output_file_name)): 
            meep.master_printf("\nWARNING: overwriting an existing HDF5 file: `%s' \n" % output_file_name)
            meep.master_printf("\tEither it remained from previous simulation, or two slices share the same definition\n")
            meep.master_printf("\tof name, space and time, in which case the export will fail. If you want to export field as a vector \n")
            meep.master_printf("\tfrom one slice, you may use e.g. `components=(meep.Ex, meep.Ey, meep.Ez)'\n\n")
        meep.all_wait()
        self.openfile = meep.prepareHDF5File(output_file_name)

    def poll(self, now):
        """ Check whether the time has come to save a slice

        The first condition is for a time span, the second one is for single snapshot.
        """
        if (((now > self.at_t[0]) and (now < self.at_t[1]) and (now-self.last_slice_time > self.min_timestep)) or 
                (now > self.at_t[0] and self.images_number==0)):
            self.images_number += 1 
            for component in self.components:
                self.field.output_hdf5(component, self.volume, self.openfile, 1) 
            self.last_slice_time = now

    def finalize(self, forcesave=True):
        if forcesave:
            self.images_number += 1 
            for component in self.components:
                self.field.output_hdf5(component, self.volume, self.openfile, 1) 
        del(self.openfile)          ## all processes must release the HDF5 file
        if meep.my_rank() == 0:        ## but postprocessing is to be done by a single process
            if self.outputgif or self.outputpng:
                meep.master_printf("Generating %d images\n" % self.images_number)
                #run_bash("cd %s; h5topng -t 0:%d -R -Zc dkbluered -a yarg %s.h5 -S 1 -o %s.png" %   ## XXX
                run_bash("cd %s; h5topng -t 0:%d -R -Zc dkbluered -a yarg %s.h5 -S 1 -o %s.png" % 
                        (self.outputdir, self.images_number-1, self.name, self.name))
            if self.outputgif: 
                meep.master_printf("Converting %d images to gif\n" % self.images_number)
                run_bash("cd %s; convert -compress None -delay 10 %s*png %s.gif" % (self.outputdir, self.name, self.name))
            if self.outputgif and not self.outputpng: 
                run_bash("cd %s; rm %s*.png" % (self.outputdir, self.name))

            # TODO if self.outputplot: 

            if self.outputvtk: 
                export_args = ""
                for component in self.components:
                    comp_name = meep.component_name(component)
                    if comp_name not in ('eps', 'cond', 'mu'): comp_name += self.real_or_imag
                    export_args += "%s.h5:%s " % (self.name, comp_name)
                flatten_time =  "-t 0" if not self.isrange(self.at_t) else "" ## avoid flat 4th dimension (otherwise vtk fails)
                run_bash ("cd %s; h5tovtk %s %s -o %s.vtk" % 
                        (self.outputdir, export_args, flatten_time, self.name))

            if not self.outputhdf: 
                run_bash("cd %s; rm %s.h5" % (self.outputdir, self.name))
    #}}}

## Obtain and process the s-parameters of the structure 
def get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, #{{{
        frequency_domain=False, frequency=None, pad_zeros=0.0, intf=[0, np.inf], Kx=0, Ky=0, eps1=1, eps2=1):
    """ Returns the frequency, s11 (reflection) and s12 (transmission) spectra
    (works for both time- and freq-domain simulation) 

    Separate forward and backward waves in time domain                   +-----------------+           
    (both electric and magnetic fields needed for this)       -- in1 ->|MP1|--+---s12--->|MP2| -- out2 --> 
         a ... inputs; b ... outputs                                   |MP1| s11         |MP2|       
         1 ... front port, 2 ... rear port (on z-axis)        <- out1 -|MP1|<-'          |MP2| <- in2=0 --
    MP1 and MP2 are monitor planes provided                              +------structure--+           

    Allows to separate forward/backward waves even under oblique incidence
    """
    ## TODO document function parameters
    ## TODO allow omitting second monitor (-> returns s12=None)

    t, Ex1 = monitor1_Ex.get_waveforms()
    t, Hy1 = monitor1_Hy.get_waveforms()
    t, Ex2 = monitor2_Ex.get_waveforms()
    t, Hy2 = monitor2_Hy.get_waveforms()

    ## Hann-window fadeout to suppress spectral leakage
    if not frequency_domain:
        for field in (Ex1, Hy1, Ex2, Hy2):
            field[t>max(t)*.8] = field[t>max(t)*.8]*(.5 + .5*np.cos(np.pi * (t[t>max(t)*.8]/max(t)-.8)/(1-.8)))

    try:
        import matplotlib
        #matplotlib.use('Agg') ## Enable plotting even in the GNU screen session?
        from matplotlib import pyplot as plt
        plt.figure(figsize=(7,6))
        plt.plot(t, abs(Ex1), label="Ex1")
        plt.plot(t, abs(Hy1), label="Hy1")
        plt.plot(t, abs(Ex2), label="Ex2", lw=1.5)
        plt.plot(t, abs(Hy2), label="Hy2")

        plt.gca().set_ylim(ymin=1e-10)
        plt.legend(prop={'size':10}, loc='upper right')
        plt.xlabel('Time'); plt.ylabel('Field amplitudes, $|E|$, $|H|$')
        plt.yscale("log")
        plt.savefig("amplitudes_time_domain.png", bbox_inches='tight')
    except:
        print "Timedomain plot failed"


    ## Obtain the electric and magnetic fields spectra
    if frequency_domain:            ## No need for FFT in frequency domain, just copy the value
        freq = np.array([frequency])
        (Ex1f, Hy1f, Ex2f, Hy2f) = (Ex1, Hy1, Ex2, Hy2)
    else:
        if pad_zeros: 
            ## Extend the data range by zeros so that FFT is efficient; (artificially prolonging stable eff param retrieval)
            target_len = 2**np.ceil(np.log(len(Ex1)*(1+pad_zeros))/np.log(2))      ## must be power of two for efficient FFT!
            append_len = target_len - len(Ex1)
            Ex1, Hy1, Ex2, Hy2  =  map(lambda x: np.append(x, np.zeros(append_len)), (Ex1, Hy1, Ex2, Hy2))

        ## Calculate the Fourier transform of the recorded time-domain waves
        numpoints = len(Ex1)
        freq = np.arange(0., int(numpoints/2)) / numpoints / (t[1]-t[0]) # take positive frequency range only
        ## TODO use np.fft.fftfreq() instead:
        ## fftfreq(signal.size, Sample spacing.[d])	Return the Discrete Fourier Transform sample frequencies.
        ## TODO Positive frequency range should be  separated just by truncation below, as 'minf => 0'

        #fftshift(x[, axes])	Shift the zero-frequency component to the center of the spectrum.
        (Ex1f, Hy1f, Ex2f, Hy2f) = map(lambda x: np.fft.fft(np.real(x))[0:int(numpoints/2)], (Ex1, Hy1, Ex2, Hy2))
        
        ## Truncate the data ranges to allowed radiating angles, and possibly to minf<freq<maxf
        truncated = np.logical_and(np.logical_and((Ky**2+Kx**2)<((2*np.pi*freq/c)**2), freq>intf[0]), freq<intf[1])
        (Ex1f, Hy1f, Ex2f, Hy2f, freq) = map(lambda x: x[truncated], (Ex1f, Hy1f, Ex2f, Hy2f, freq))


    ## Diagnostics: plot frequency-domain data
    try:
        plt.figure(figsize=(7,6))
        plt.plot(freq, abs(Ex1f), label="Ex1")
        plt.plot(freq, abs(Hy1f), label="Hy1")
        plt.plot(freq, abs(Ex2f), label="Ex2")
        plt.plot(freq, abs(Hy2f), label="Hy2")
        plt.yscale("log");   plt.gca().set_ylim(ymin=1e-8)
        plt.yscale("log");   plt.gca().set_ylim(ymin=1e-8)
        plt.xlim(0, np.max(freq))
        plt.legend(prop={'size':10}, loc='upper right')
        plt.xlabel('Frequency'); plt.ylabel('Field amplitudes, $|E|$, $|H|$')
        plt.savefig("amplitudes_freq_domain.png", bbox_inches='tight')
    except:
        print "Raw freq-domain plot failed", sys.exc_info()[0]


    ## Prepare the angles at which the wave propagates (dependent on frequency, Kx and Ky)
    ## Separate the forward and backward wave in frequency domain 
    ##    (Efield+Hfield)/2 ->    forward wave amplitude, 
    ##    (Efield-Hfield)/2 ->    backward wave amplitude
    in1, out1 =  (Ex1f+Hy1f)/2, (Ex1f-Hy1f)/2 ## old: works only for perp. incidence beta0=0
    in2, out2 =  (Ex2f-Hy2f)/2, (Ex2f+Hy2f)/2
    #beta0 = np.arcsin((Kx**2+Ky**2)**.5 / (2*np.pi*freq/c))
    #in1, out1 =  (Ex1f+Hy1f/np.cos(beta0))/2, (Ex1f-Hy1f/np.cos(beta0))/2 ## old: works only for monitors placed in vacuum
    #in2, out2 =  (Ex2f-Hy2f/np.cos(beta0))/2, (Ex2f+Hy2f/np.cos(beta0))/2


    # **********************
    # amplitude squared is double Poynting vector (in any medium):   2 S = E H = E**2/Z = E**2 * eps**.5 / mu**.5 = E**2 * c * eps

    # if there are no losses, we assert:   |in1|**2 + |in2|**2 = |out1|**2 + |out2|**2
    # -> which trivially leads to:    r**2 + t**2 = 1   (if in2=0)

    # **********************

    #n1, n2 = eps1**.5,  eps2**.5
    #z1, z2 = eps1**-.5, eps2**-.5
    #Ktot1 = (2*np.pi*freq*n1/c)
    #angles1 = np.arcsin((Kx**2+Ky**2)**.5 / Ktot1)
    #in1 =  (Ex1f/z1  +Hy1f/np.cos(angles1))/2/eps1**.5 * n1**.5 
    #out1 = (Ex1f/z1  -Hy1f/np.cos(angles1))/2/eps1**.5 * n1**.5
    #
    #Ktot2 = (2*np.pi*freq*n1/c)
    #angles2 = np.arcsin((Kx**2+Ky**2)**.5 / Ktot2)
    #in2  = (Ex2f/z2  -Hy2f/np.cos(angles2))/2/eps2**.5 * n2**.5
    #out2 = (Ex2f/z2  +Hy2f/np.cos(angles2))/2/eps2**.5 * n2**.5
    
    #in2, out2 =  (Ex2f      -Hy2f/np.cos(beta0mon2))/2, (Ex2f  +Hy2f/np.cos(beta0mon2))/2
    ## Todo optimize cos(arcsin x ) = sqrt(1-x**2)

    ## Diagnostics: Plot spectral profile
    try:
        plt.figure(figsize=(7,6))
        plt.plot(freq, abs(in1), label="in1")
        plt.plot(freq, abs(out1), label="out1")
        plt.plot(freq, abs(in2), label="in2")
        plt.plot(freq, abs(out2), label="out2")
        plt.xlim(0, np.max(freq))
        plt.legend(prop={'size':10}, loc='lower left')
        plt.xlabel('Frequency'); plt.ylabel('Transmitted amplitude')
        #plt.title('Frequency-domain wave amplitudes')
        plt.yscale("log")
        plt.savefig("amplitudes_spectra.png", bbox_inches='tight')
    except:
        print "Wave amplitude freq-domain plot failed", sys.exc_info()[0]

    ## Get the s-parameters 
    s11 = out1 / in1
    s12 = out2 / in1

    ## Return the S-parameters (i. e. complex reflection and transmission)
    return freq, s11, s12, "#x-column freq\n#column |r|\n#column r phase\n#column |t|\n#column t phase\n"
#}}}
def get_phase(complex_data):#{{{
    """ Unwraps and shifts the phase from Fourier transformation """
    if len(complex_data) <= 1: return np.angle(complex_data)
    phase = np.unwrap(np.angle(complex_data))
    center_phase = phase[min(5, len(phase)-1)] ## 5 is chosen to avoid zero freq.
    return phase-(round(center_phase/2/np.pi)*2*np.pi)
#}}}
class AmplitudeMonitorPlane():#{{{
    """ Calculates an average of electric field and perpendicular magnetic field.

    I asked for a similar field-averaging function built in MEEP, but it seems not to be implemented yet.
    http://www.mail-archive.com/meep-discuss@ab-initio.mit.edu/msg04447.html
    A field averaging routine written in Python is used, which is very inefficient. 

    Note this implementation requires the planes are in vacuum (where impedance = 1.0)
    """
    def __init__(self, comp=None, size_x=None, size_y=None, z_position=None, Kx=0, Ky=0):
        self.comp=comp
        self.size_x = size_x
        self.size_y = size_y
        self.z_position = z_position
        self.Kx = Kx
        self.Ky = Ky

        self.t = []
        self.waveform = []

    def average_field(self, field):
        """
        Average field component in some plane, return amplitudes 
        This function is ineffective - it should be implemented in C++ in meep itself

        5x5 grid is usually optimal (no visible difference between 10x10 grid and 5x5 grid)

        TODO:  This class implements a workaround for unavailable amplitude averaging in python-meep.
        This implementation is ineffective and inflexible, but one would have to edit the MEEP source otherwise.
        """
        xcount, ycount = (3, 3)
        field_sum = 0 
        # The mode function has the form of an oblique plane wave
        #for x in [x0*self.size_x/xcount+(self.size_x/2/xcount)-self.size_x/2 for x0 in range(xcount)]:
            #for y in [y0*self.size_y/ycount+(self.size_y/2/ycount)-self.size_y/2 for y0 in range(ycount)]:
                #field_sum += (field.get_field(self.comp, meep.vec(x, y, self.z_position)) *  
                                #np.exp(1j*(self.Kx*x + self.Ky*y)) )
        #return field_sum/(xcount*ycount)

        ## New way (removes explicit cycle, few percent faster)
        xr = [x0*self.size_x/xcount+(self.size_x/2/xcount)-self.size_x/2 for x0 in range(xcount)]
        yr = [y0*self.size_y/ycount+(self.size_y/2/ycount)-self.size_y/2 for y0 in range(ycount)]
        xm, ym = np.meshgrid(xr,yr)
        points = zip(xm.flatten(), ym.flatten())
        sum_ = sum(map(lambda pos: field.get_field(self.comp, meep.vec(pos[0], pos[1], self.z_position)), points))
        return sum_/(xcount*ycount)

        ## Yet newer way
        #xr = [x0*self.size_x/xcount+(self.size_x/2/xcount)-self.size_x/2 for x0 in range(xcount)]
        #yr = [x0*self.size_x/xcount+(self.size_x/2/xcount)-self.size_x/2 for x0 in range(xcount)]
        #xm, ym = np.meshgrid(xr,yr)
        #v = meep.vec(0,0, self.z_position)
        #points = zip(xm.flatten(), ym.flatten())
        #sum_ = sum(map(lambda pos: _meep_mpi.fields_get_field(field, self.comp, v), points))
        #sum_ = sum(map(lambda pos: _meep_mpi.fields_get_field(field, self.comp, _meep_mpi.new_vec(0,0,0)), points))
        #cp =self.comp
        #zp = self.z_position
        #sum_ = sum(map(lambda pos: _meep_mpi.fields_get_field(field, cp, _meep_mpi.new_vec(pos[0], pos[1], zp)), points))
        #return sum_/(xcount*ycount)
        


    #def NEW_average_field_xy_plane(field, component, zpos, model): ## TODO use the internal integration of fields by MEEP
        # TODO rewrite:
        # (fields &f, linear_integrand_data &d, const volume &v, component cgrid)
        # f.integrate(0, 0, linear_integrand, (void *) &d, v)
        #integrate(meep::fields *,int,meep::component const *,meep::field_function,void *,meep::volume const &,double *)
        #v = meep.volume(
                #meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness), 
                #meep.vec(model.size_x/2, model.size_y/2, -model.size_z/2+model.pml_thickness))
        #return field.integrate(1, [component], meep.one_vec, [], v)
      #Possible C/C++ prototypes are:
        #meep::fields::integrate(int,meep::component const *,meep::field_function,void *,meep::volume const &,double *)
        #meep::fields::integrate(int,meep::component const *,meep::field_function,void *,meep::volume const &)
        #meep::fields::integrate(int,meep::component const *,meep::field_rfunction,void *,meep::volume const &,double *)
        #meep::fields::integrate(int,meep::component const *,meep::field_rfunction,void *,meep::volume const &)
        #fields::integrate(int num_fvals, const component *components,
        #field_function integrand,
        #void *integrand_data_,
        #const volume &where,
        #double *maxabs
    
    def record(self, field=None):
        """ 
        Useful for time-domain simulation only
        """
        self.t.append(field.time()/c)
        self.waveform.append(self.average_field(field))

    def get_waveforms(self):
        """ Return the recorded waveform (for time domain simulation only) """
        if len(self.t) <= 1:
            t, result_wform = np.array(self.t), np.array(self.waveform)
        else:
            t = np.array(self.t[:-1])
            ## The FDTD calculation introduces half-step time shift between Ex and Hy. Compensated by averaging the Hy field
            ## with its value in a next timestep. The error is reduced from O1 to O2.
            ## See http://ab-initio.mit.edu/wiki/index.php/Synchronizing_the_magnetic_and_electric_fields
            if meep.is_magnetic(self.comp) or meep.is_B(self.comp):
                result_wform = np.array(self.waveform[:-1])/2. + np.array(self.waveform[1:])/2.
            else: 
                result_wform = np.array(self.waveform[:-1])
            
        return t, result_wform 
         ## time, 
        ## TODO this will have to be modified in order to account for oblique incidence
        ## TODO take into account the medium impedance (... divide the Hfield)
#}}}


## === Experimental zone ===
""" TODOs:#{{{
    * replace the classes of AmplitudeMonitorPlane and AmplitudeMonitorPoint 
        with a class AmplitudeMonitor (accepts shape={meep.vec, meep.volume}
        in __init__, it defines a list of one or more points the field is summed over
        these points contain 
            1) the position (meep.vec object), 
            2) amplitude as a complex number, this allows to define a plane that detects oblique incidence,
            3) meep component
        so a port recording a forward wave with +amplitude and backward wave with -amplitude
    *  ... or is it useless and everything can be extracted from a Slice? (It would avoid slow get_field cycles, 
        and would make the scripts a bit more compatible with B-Calm.)

       Slice.finalize() 
                returns the data as numpy.array (use `pytables' to load the HDF5 file back)
                if the Slice.component is magnetic, interpolate two adjacent numbers
                if the Slice.component is electric, just remove first number to get the same time axis

       Slice.data()     gives this numpy.array
       simpleFFT(Slice.data())  computes FFT along the time axis (this is cool, because FFT 
                 should inherently return the same value in freq-domain simulation)
                 accepts a parameter of whether to pad time axis to 2**n; true by default

       Slice.axes() gives a dictionary as such  {'x':np.array, 'y':np.array, 't':np.array}
       Slice.axesFFT()          dtto            {'x':np.array, 'y':np.array, 'f':np.array}

    *  SpectralModeAmplitude(simpleFFT(Slice.data()), ModeFunction, kwargs_for_ModeFunction)  -->  gives amplitude in the mode
            Modefunction(r, freq, **kwargs_for_ModeFunction)  may be dependent on frequency to accomodate eg. 
                1) for oblique plane wave (combines Hy*cos + Hz*sin)
                2) for modes in dielectric waveguides etc. which change their profile with frequency a bit
            returns forward and backward amplitudes
    *  Slices_to_Amplitude(k_vector, slices)
        k_vector  .. (kx,ky,kz)
        slices     ... {"Ex": np.array, "Ey": ...    "Hz"}
                (generally the function will accept all 6 components, 

        In my case of metamaterial research, a snippet of s-parameter computation will be like this:
            k1 = (model.Kx, model.Ky, ((2*np.pi*freq*slice_refr_index)**2 - Kx**2 - Ky**2)**.5)
            k2 = (model.Kx, model.Ky, ((2*np.pi*freq*slice_refr_index)**2 - Kx**2 - Ky**2)**.5)
            amplitude_front_in, amplitude_front_out = Slices_to_Amplitude(k1, SliceEx, ...
            amplitude_rear_in,  amplitude_rear_out  = ...
            s11 = amplitude_front_out / amplitude_front_in
            s12 = amplitude_rear_out / amplitude_front_in
"""#}}}
def lorentzian_unstable_check_new(model, dt, quit_on_warning=True): #{{{ ## disused
    """
    Experimental routine
    """
    for mat in model.materials:
        eps_ts = analytic_eps(mat, 1/dt/np.pi)  
        if np.real(eps_ts)<0: ## <--- this is WRONG, corrected newly in test_materials()
            meep.master_printf("Warning: for material '%s', the permittivity is negative at the critical frequency eps(1/pi/dt)=eps(%g)=%s.\n" % \
                        (mat.name, 1/dt/np.pi, eps_ts.__str__()));
            if quit_on_warning: quit()
        for pol in mat.pol:
            omega_0, gamma = pol['omega'], pol['gamma']   ## (non-angular!) frequency of the oscillator and damping rate
            if (omega_0 > gamma/2):
                z2 = np.sqrt(gamma*gamma + 4*omega_0*omega_0)/2
            else:
                z2 = gamma/2 + np.sqrt(gamma*gamma - 4*omega_0*omega_0)/2
            if z2 > 1/dt/np.pi:     ## <----- this is QUESTIONABLE, todo: check the actual stability when gamma is increased 
                meep.master_printf("Warning: for material '%s', the oscillator pole having magnitude |z|=%g will be probably unstable when 1/pi/dt=%g.\n" % \
                        (mat.name, z2, 1/dt/np.pi));
                if quit_on_warning: quit()
#}}}

def save_s_params_old(name="output.dat", freq=None, s11=None, s12=None, model=None, polar_notation=True, truncate=True): #{{{     
    # TODO rewrite simulations to use save_txt instead, and delete this mercilessly
    """ Saves the s-parameters to a file, including comments """
    with open(model.simulation_name+".dat", "w") as outfile: 
        outfile.write(model.parameterstring)

        ## Truncate the data sets
        if truncate:
            mask = np.logical_and(freq>interest_freq[0], freq<interest_freq[1])
            (s11, s12, freq) = map(lambda arr: arr[mask], (s11, s12, freq))

        if polar_notation:
            ## Convert to polar notation
            s11amp, s12amp, s11phase, s12phase = abs(s11), abs(s12), get_phase(s11), get_phase(s12)
            ## Save polar
            outfile.write("#x-column Frequency [Hz]\n#Column Reflection amplitude\n#Column Reflection phase\n" + \
                        "#Column Transmission amplitude\n#Column Transmission phase\n")
            np.savetxt(outfile, zip(freq, s11amp, s11phase, s12amp, s12phase), fmt="%.8e")
        else:
            ## Save cartesian
            # TODO should save in such format that PKGraph understands complex values
            outfile.write("#x-column Frequency [Hz]\n#Column Reflection Re\n#Column Reflection Im\n" + \
                        "#Column Transmission Re\n#Column Transmission Im\n")
            np.savetxt(outfile, zip(freq, s11.real, s11.imag, s12.real, s12.imag), fmt="%.8e")
#}}}
def load_rt_old(filename, layer_thickness=None, plot_freq_min=None, plot_freq_max=None, truncate=True): #{{{   # TODO rename this to load_s_params
    """ Loads the reflection and transmission spectra and simulation settings 

    Returns:
    * frequency axis
    * reflection s11 and transmission s12 as complex np arrays

    Compatible with the program used in our laboratory (PKGraph), uses polar notation for complex data: 
    * parameters in header like: #param name,value
    * column identification like: #column Ydata
    * data columns in ascii separated by space
    Expects polar data with columns: frequency, s11 ampli, s11 phase, s12 ampli, s12 phase
    """
    with open(filename+'.dat') as datafile:
        for line in datafile:
            if line[0:1] in "0123456789": break         # end of file header
            value = line.replace(",", " ").split()[-1]  # the value of the parameter will be separated by space or comma
            if ("layer_thickness" in line) and (layer_thickness == None): d = float(value)
            if ("plot_freq_min" in line) and (plot_freq_min == None): plot_freq_min = float(value)
            if ("plot_freq_max" in line) and (plot_freq_max == None): plot_freq_max = float(value)
    xlim = (plot_freq_min, plot_freq_max)
    (freq, s11amp, s11phase, s12amp, s12phase) = \
            map(lambda a: np.array(a, ndmin=1), np.loadtxt(filename+".dat", unpack=True)) 

    ## Limit the frequency range to what will be plotted (recommended)
    #TODO better: 
    #truncated = np.logical_and(freq>minf, freq<maxf) 
    #(a, b, c, freq) = map(lambda x: x[truncated], (a, b, c, freq))
    if truncate:
        (d0,d1) = np.interp((plot_freq_min, plot_freq_max), freq, range(len(freq)))
        (freq, s11amp, s11phase, s12amp, s12phase) = \
                map(lambda a: a[int(d0):int(d1)], (freq, s11amp, s11phase, s12amp, s12phase))
    return freq, s11amp, s11phase, s12amp, s12phase, xlim, (d, plot_freq_min, plot_freq_min)
#}}}

class AmplitudeMonitorPoint(AmplitudeMonitorPlane):#{{{
    """ Calculates an average of electric field and perpendicular magnetic field.
    """
    def __init__(self, Ecomp=None, Hcomp=None, pos=None):
        self.Ecomp=Ecomp
        self.Hcomp=Hcomp
        self.pos = pos          ## where to record the field
        self.t = []
        self.Efield = []
        self.Hfield = []

    def get_amplitude(self, field, component):
        """ Record field in some point. No averaging here, but (inherited) time recording is pretty useful. """
        return field.get_field(component, self.pos)
#}}}
        
"""
        
== Possible errors and clues ==
    HDF5-DIAG: Error detected in HDF5 (1.8.4-patch1) MPI-process 0:
      #000: ../../../src/H5F.c line 1514 in H5Fopen(): unable to open file
    ---> perhaps you defined two slices the same?

mpirun has exited due to process rank * with PID **** on *
    ---> this is harmless

"""
