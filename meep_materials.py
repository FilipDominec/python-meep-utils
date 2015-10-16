#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
meep_materials.py
This script contains definitions of materials suitable for the FDTD algorithm.
However, only the materials from the sections "Generic" and "Prepared for FDTD" are
ready to be fed to MEEP, and they are valid only for a given spectral range. 

On the contrary, the materials in the "Realistic" section are defined for whole
spectrum from microwave to optical range, so the number of Lorentzian oscillators
would have to be manually reduced to make the simulation run effectively.

About this script:
 * Written in 2012-2015 by Filip Dominec (dominecf at the server of fzu.cz)
 * Being distributed under the GPL license, this script is free as speech after five beers. 
 * You are encouraged to use and modify it as you need. Feel free to write me if needed.
 * Hereby I thank to the MEEP/python_meep authors and people of meep mailing list who helped me a lot.
 
See also:
    http://f.dominec.eu/meep/              -- my MEEP page
    http://f.dominec.eu/eps/               -- permittivity spectra plot into graphs

TODO:
    * implement material self-adjustment for stability 
       1) wipe out all oscillators above f_c (sum Δε into ε_infty)
       2) adjust the Drude term to be stable when f_c < omega_p/2π
       3) re-express all oscillators above FRoI as one oscillator (is it possible? what are the summing rules
                for lorentzians with different frequency and losses?)
       4) do something with that Debye term can be unstable
    * make a manifest of materials and their source files; 
    * add the plot_eps.py file and its directory into the project??

Last edited: 2015-02-19
"""

import numpy as np
from scipy.constants import epsilon_0, c, pi
percm = c/1e-2



## -- Generic --
class material_dielectric():#{{{
    """
    This model defines a simple nondispersive dielectric. 
    The predefined permittivity of 2 holds e. g. for polyethylene, teflon etc. and other materials that 
    have a flat index of refraction ~ 1.4 both at microwave and optical frequencies
    """
    def __init__(self, where=None, eps=2., loss=0, chi2=0, chi3=0):
        self.eps = eps*(1-loss)
        ## TODO: make the loss effect well defined (for some frequency), and write 
        if loss != 0:
            self.pol = [
                    {'omega':5e12, 'gamma':10e12,  'sigma':100*loss},            ## TODO make clearer
                    ]
        else:
            self.pol = []
        self.name = "Lossy dielectric"
        #self.chi2 = chi2
        #self.chi3 = chi3
        self.where = where
#}}}
# TODO def MakeDrudian() 
class material_DrudeMetal():#{{{
    """ Defines a generic metal with a Drude model

    In Drude model, the electrons are expected to move freely with some inertia, which
    together with their density determines the plasma frequency omega_p.
    Losses are introduced by letting them to randomly scatter with the frequency gamma.
    The relative permittivity eps_r(omega) of metal can then be calculated as:

        eps_r(omega) = 1 - omega_p**2 / (omega**2 - i*omega*gamma), 
        
    Low-frequency limit of permittivity (eps = eps' + 1j * eps''):
            eps'   --->     - omega_p**2/gamma**2,  
            eps''  --->     omega_p / (omega * gamma),     
    and the high-frequency limit of permittivity:
            eps'   --->     1 - omega**2/omega_p**2         
            eps''  --->     1/omega**3                      
    """

    #def __init__(self, where=None, lfconductivity=15e6, f_c=1e14): XXX
    def __init__(self, where=None, lfconductivity=15e3, f_c=1e14, gamma_factor = .5, epsplus=0):
        """
        At microwave frequencies (where omega << gamma), the most interesting property of the metal 
        is its conductivity sigma(omega), which may be approximated by a nearly constant real value:
            
            conductivity = epsilon/1j * frequency * eps0 * 2pi

        Therefore, for convenience, we allow the user to specify the low-frequency conductivity sigma0 and compute
        the scattering frequency as:

        Values similar to those of aluminium are used by default:
            dielectric part of permittivity: 1.0, 
            plasma frequency = 3640 THz
            low-frequency conductivity: 40e6 [S/m]

        The gamma_factor is required for numerical stability. It has always to be lower than the timestepping frequency f_c, 
        but sometimes the simulation is unstable even though and then it has to be reduced e.g. to 1e-2 or to even smaller value.
        """

        ## Design a new metallic model. We use 'omega_p' and 'gamma' in angular units, as is common in textbooks
        ## We need to put the scattering frequency below fc, so that Re(eps) is not constant at f_c (it would hinder
        ## our trick that ensures FDTD stability)
        self.gamma      = f_c * gamma_factor * 2*np.pi
        ## TODO: test when unstable -> reduce self.gamma by 100 -> test again -> finally write automatic selection of value

        ## Knowing the scattering frequency gamma, the virtual plasma frequency is now determined by lfconductivity
        omega_p = np.sqrt(self.gamma * lfconductivity / epsilon_0)

        ## The following step is required for FDTD stability: 
        ## Add such an high-frequency-epsilon value that shifts the permittivity to be positive at f_c
        ## If self.eps>1, the frequency where permittivity goes positive will generally be less than f_p
        self.eps = max((omega_p/2/np.pi / f_c)**2,   1.) + epsplus
        #print 'self.eps', self.eps
        #print "gamma =", self.gamma
        #print 'omega_p = ', omega_p
        #print 'LFC = %.3e' % (omega_p**2 * epsilon_0 / self.gamma)

        ## Feed MEEP with a Lorentz oscillator of arbitrarily low frequency f_0 so that it behaves as the Drude model
        omega_0 = .1           
        self.pol = [
                {'omega': omega_0/(2*np.pi), 'gamma': self.gamma/(2*np.pi), 'sigma': (omega_p/omega_0)**2}, # (Lorentz) model
                ## Note: meep also uses keywords 'omega' and 'gamma', but they are non-angular units
                ]
        self.name = "Drude metal for up to %.2g Hz" % f_c
        self.where = where
#}}}

## -- Prepared for FDTD simulations in the terahertz range (obsoleted) --
class material_Metal_THz():#{{{ ## Obsoleted
    """
    This model, by default, defines a metal similar to aluminium. Its parameters are roughly: 
        dielectric eps=1., 
        plasma frequency = 3640 THz [see e. g. Pendry1999] and 
        losses gamma=25e12 (roughly by electron scattering time)

    However, for simulations e. g. at 1 THz, this gives complex epsilon of (-21000 + 500000j), which may be correct,
    but it makes the simulation unstable. (This can be detected soon, as computation with gigantic numbers or NaNs 
    gets few times slower.)
    
    We therefore may use a "pseudometal" with reduced plasma frequency, which has to be _lower_ than the frequency of
    timestepping. This still gives correct results for simulations in terahertz range. 

    Quite a good metal model (up to 2 THz) is: 'omega':1.39e+9, 'gamma':1.0e+10, sigma=5e8
    For sigma over approx 1e10 (depending also on resolution) the simulation diverges.
    """
    def __init__(self, where, sigmafactor=1):
        self.eps = 50. 
        omega0 = 1e-5           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        omega_p = 3640e12       ## plasma frequency of aluminium
        #sigmafactor  = 30e12**2 /omega_p**2
        #sigmafactor  = 30e12**2 /omega_p**2
        self.pol = [
                #{'omega':1.39e9, 'gamma': 1.0e10, 'sigma':2.6e9}, 
                #{'omega': omega0, 'gamma': 25e12, 'sigma': sigmafactor * omega_p**2 / omega0**2}, # (Lorentz) model for Al
                #{'omega': omega0, 'gamma': 1e12, 'sigma': sigmafactor * omega_p**2 / omega0**2}, # stable? model at THz

                {'omega':1.39e9, 'gamma': 1e10, 'sigma':1e10}, # Original value, use sigmafactor=1
                #{'omega':1.39e9, 'gamma': 1e13, 'sigma':1e10}, # Lossy
                ]
        self.name = "Metal with reduced carrier density (for stable low-res simulations in THz range)"
        self.where = where
#}}}
class material_Metal_THz_resistive():#{{{ ## Obsoleted
    """
    This model, by default, defines a metal similar to aluminium. Its parameters are roughly: 
        dielectric eps=1., 
        plasma frequency = 3640 THz [see e. g. Pendry1999] and 
        losses gamma=25e12 (roughly by electron scattering time)

    However, for simulations e. g. at 1 THz, this gives complex epsilon of (-21000 + 500000j), which may be correct,
    but it makes the simulation unstable. (This can be detected soon, as computation with gigantic numbers or NaNs 
    gets few times slower.)
    
    We therefore may use a "pseudometal" with reduced plasma frequency, which has to be _lower_ than the frequency of
    timestepping. This still gives correct results for simulations in terahertz range. 

    Quite a good metal model (up to 2 THz) is: 'omega':1.39e+9, 'gamma':1.0e+10, sigma=5e8
    For sigma over approx 1e10 (depending also on resolution) the simulation diverges.
    """
    def __init__(self, where, sigmafactor=1):
        self.eps = 50. 
        omega0 = 1e-5           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        omega_p = 3640e12       ## plasma frequency of aluminium
        #sigmafactor  = 30e12**2 /omega_p**2
        #sigmafactor  = 30e12**2 /omega_p**2
        self.pol = [
                #{'omega':1.39e9, 'gamma': 1.0e10, 'sigma':2.6e9}, 
                #{'omega': omega0, 'gamma': 25e12, 'sigma': sigmafactor * omega_p**2 / omega0**2}, # (Lorentz) model for Al
                #{'omega': omega0, 'gamma': 1e12, 'sigma': sigmafactor * omega_p**2 / omega0**2}, # stable? model at THz

                {'omega':1.39e9, 'gamma': 1e12, 'sigma':6e9}, # Original value, use sigmafactor=1
                #{'omega':1.39e9, 'gamma': 1e13, 'sigma':1e10}, # Lossy
                ]
        self.name = "Metal with reduced carrier density (for stable low-res simulations in THz range)"
        self.where = where
#}}}
class material_TiO2_THz():#{{{  ## for THz application only
    """
    This model defines the TiO2 material for 0 - 2 THz frequencies as a lossy dielectric. 
    For instance, at 500 GHz, epsilon = 94 + 1.2j
    The epsilon=22. constant is the remaining high-frequency part of permittivity (from [Baumard1977]).
    Polarizability resonances to obtain absorption [Baumard1977] with following changes:
        1) We used only one resonance at 5.6 THz, as the other should be negligible
        2) The resonance width was reported to vary from 8e11 to 1.6e12 Hz according to composition [Baumard1977].
            We used 1.05e12 to match our experimental data for losses in bulk TiO2.

    """
    def __init__(self, where=None):
        self.eps = 12.
        self.pol = [
                #{'omega':5.67e12, 'gamma':1.05e12, 'sigma':70.},            ## strongest optical phonon resonance in THz range
                {'omega':5.67e12, 'gamma':1.05e12, 'sigma':50+30.},            ## strongest optical phonon resonance in THz range
                ]
        self.name = "polycrystalline TiO2 (rutile)"
        self.where = where
#}}}
class material_TiO2_THz_HIEPS():#{{{  ## for THz application only
    """
    This model defines the TiO2 material for 0 - 2 THz frequencies as a lossy dielectric. 
    For instance, at 500 GHz, epsilon = 94 + 1.2j
    The epsilon=22. constant is the remaining high-frequency part of permittivity (from [Baumard1977]).
    Polarizability resonances to obtain absorption [Baumard1977] with following changes:
        1) We used only one resonance at 5.6 THz, as the other should be negligible
        2) The resonance width was reported to vary from 8e11 to 1.6e12 Hz according to composition [Baumard1977].
            We used 1.05e12 to match our experimental data for losses in bulk TiO2.
    """
    def __init__(self, where=None):
        self.eps = 22.*2
        self.pol = [
                #{'omega':5.67e12, 'gamma':1.05e12, 'sigma':70.},            ## strongest optical phonon resonance in THz range
                {'omega':5.67e12, 'gamma':1.05e12, 'sigma':70.*2},            ## strongest optical phonon resonance in THz range
                ]
        self.name = "polycrystalline TiO2 (rutile)"
        self.where = where
#}}}
class material_TiO2_NC_THz():#{{{  ## for THz application only
    """
    TiO2 used in Navarro-Cia publication:

    "The TiO2 particle is modelled by a perfect 
    sphere of radius r = 15 um and electric
    permittivity r = 94 + 2.35i (tan = 0.025)"
    """
    def __init__(self, where=None):
        self.eps = 21.8
        self.pol = [
                #{'omega':5.67e12, 'gamma':1.05e12, 'sigma':70.},            ## strongest optical phonon resonance in THz range
                {'omega':5.67e12, 'gamma':1.05e12*(2.35/2.43), 'sigma':70.},            ## strongest optical phonon resonance in THz range
                ]
        self.name = "polycrystalline TiO2 (rutile)"
        self.where = where
#}}}
class material_STO_THz():#{{{  ## for THz application only
    """
    STO

    1st phonon shifts with temperature:
    4K: 250 GHz, 22000
    90K: 850 GHz, 1400
    300K: 2300 GHz, 310
    """
    def __init__(self, where=None):
        self.eps = 100.
        self.pol = [
                #{'omega':2e12, 'gamma':.7e12, 'sigma':210.},
                {'omega':2e12, 'gamma':.3e12, 'sigma':210.}, ## XXX 
                ]
        self.name = "Strontium titanate (T = 300 K) for THz"
        self.where = where
#}}}
class material_STO_hiloss():#{{{  ## for THz application only
    """
    STO

    1st phonon shifts with temperature:
    4K: 250 GHz, 22000
    90K: 850 GHz, 1400
    300K: 2300 GHz, 310
    """
    def __init__(self, where=None):
        self.eps = 100.
        self.pol = [
                {'omega':2e12, 'gamma':1e12, 'sigma':210.},
                ]
        self.name = "Strontium titanate (T = 300 K) for THz, double loss"
        self.where = where
#}}}
class material_STO_HTChen():#{{{
    """
    This model defines strontium titanate that fits the HT Chen  structure
    """
    def __init__(self, where=None, eps=2., loss=.5):
        self.eps = eps*(1-loss)
        ## TODO: make the loss effect well defined (for some frequency)
        if loss != 0:
            self.pol = [
                    {'omega':5.67e12, 'gamma':4.05e12, 'sigma':eps*loss},            ## XXX test
                    {'omega':2e12, 'gamma':2e13, 'sigma':100*loss},            ## XXX test
                    {'omega':6e12, 'gamma':6e13, 'sigma':100*loss},            ## XXX test
                    ]
        else:
            self.pol = []
        self.name = "Common dielectric (PE, PTFE etc.)"
        self.where = where
#}}}
#class material_Sapphire_THz():#{{{  ## for THz application only
#}}}


## -- Realistic (simulations may require to make a copy and optimize) --
class material_STO():#{{{  
    """
    Strontium titanate

    1st phonon shifts with temperature
    4K: 250 GHz, 22000
    90K: 850 GHz, 1400
    300K: 2300 GHz, 310
    """
    def __init__(self, where=None):
        #self.eps = 100.
        self.eps = 1.
        self.pol = [
                {'omega':2.3e12, 'gamma':.4e12*1.000, 'sigma':310},
                {'omega':15.5e12, 'gamma':5e12, 'sigma':6.5},
                {'omega':1e15, 'gamma':1e13*1.000, 'sigma':4},
                ]
        self.name = "Strontium titanate (T = 300 K)"
        self.where = where
#}}}
class material_Sapphire():#{{{
    """

    Drude-Lorentz model (to be checked with references below)
    Thomas M. E. et al.: "Frequency and temperature dependence of the refractive index of sapphire", 
    Infrared Physics & Technology, Volume 39, Issue 4, 1998
    (from http://www.sciencedirect.com/science/article/pii/S1350449598000103)

    Compare with 
    Schubert M. et al.: "Infrared dielectric anisotropy and phonon modes of sapphire"
    Phys. Rev. B 61, 8187–8201 (2000)



    http://www.mt-berlin.com/frames_cryst/descriptions/sapphire.htm:
    Refractive index at 0.532 µm	n0=1.7717, ne=1.76355

    Schott PDF:
    Dielectric Constant 11.5 (103 to 109 Hz, 25 °C) parallel to c-axis (extraordinary)
                          9.3 (103 to 109 Hz, 25 °C) perpendicular to c-axis (ordinary)
    Loss Tangent 8.6 x 10–5 (@1010 Hz, 25 °C) parallel to c-axis
                     3.0 x 10–5 (@1010 Hz, 25 °C) perpendicular to c-axis

    http://www.tydexoptics.com/products/thz_optics/thz_materials/:
    absorption spectrum to be fit
    
    """
        
    def __init__(self, where=None, ordinary=1):
        self.eps = 1.
        extraordinary = 1 - ordinary
        self.pol = [
                #Ordinary (tg delta ~ 4e-6 @ 1e10 Hz)
                {'omega':1.1517E+13, 'gamma':1.2669E+11, 'sigma':3.3000E-01*ordinary},
                #{'omega':1.3154E+13, 'gamma':7.8923E+10, 'sigma':2.7880E+00*ordinary},
                {'omega':1.3154E+13, 'gamma':.5E+12      , 'sigma':2.7880E+00*ordinary}, ## increases losses in THz
                #{'omega':1.7029E+13, 'gamma':2.0435E+11, 'sigma':2.9800E+00*ordinary},
                {'omega':1.7029E+13, 'gamma':.5E+12     , 'sigma':2.9800E+00*ordinary},
                {'omega':1.8989E+13, 'gamma':1.8989E+11, 'sigma':1.4500E-01*ordinary},
                {'omega':2.4264E+13, 'gamma':3.8094E+12, 'sigma':1.8500E-02*ordinary},
                {'omega':2.5691E+15, 'gamma':2.5691E+10, 'sigma':6.5069E-01*ordinary},
                {'omega':4.1245E+15, 'gamma':4.1245E+10, 'sigma':1.4314E+00*ordinary},

                #Extraordinary (tg delta ~ 1.1e-5 @ 1e10 Hz)
                {'omega':1.1973E+13, 'gamma':2.3946E+11, 'sigma':6.7850E+00*extraordinary},
                {'omega':1.7502E+13, 'gamma':6.1259E+11, 'sigma':1.5600E+00*extraordinary},
                {'omega':1.9600E+13, 'gamma':1.1760E+12, 'sigma':1.6000E-02*extraordinary},
                {'omega':2.4636E+15, 'gamma':2.4636E+10, 'sigma':5.5069E-01*extraordinary},
                {'omega':4.0484E+15, 'gamma':4.0484E+10, 'sigma':1.5042E+00*extraordinary},
                ]
        self.name = "Sapphire"
        self.where = where
#}}}
class material_TiO2():#{{{
    """
    This model defines the rutile (titanium dioxide, TiO2) for 0 - 500 THz frequencies as a lossy dielectric. 

    Rutile is anisotropic. By default a polycrystalline averaging is assumed. Set the parameter:
        extraordinary=1.0 for extraordinary ray
        extraordinary=0.0 for ordinary ray

    If you run simulations in the 0-4 THz range, feel free to substitute the second and higher resonances by 
        eps=7 or something.
        The magnitude of the first resonance is more importand, one has to fine tune it.
    Additionally, the width of the resonance at omega=5.6 THz can be changed from gamma=810 GHz to higher values 
        to account for higher losses due to impurities etc. (For example, 1.05e12 better matched THz experimental data.)
    The ultraviolet absorption could not be exactly modelled using finite number of lorentzian oscillators. Two 
        of them are used, but note they introduce unrealistic high damping in the optical region.
    """
    def __init__(self, where=None, extraordinary=.33):
        self.eps = 1
        self.pol = [
                # Optical phonon resonance in THz range from Baumard1977, manually tuned to better fit experiment
                {'omega':5.67e12, 'gamma':.85e12, 'sigma':50+90*extraordinary},      
                {'omega':11.43e12, 'gamma':.495e12, 'sigma':0.9*extraordinary},            
                {'omega':15.24e12, 'gamma':.72e12, 'sigma':2.6*extraordinary},            
                {'omega':1.0e15, 'gamma':.15e15, 'sigma':2.8 + .5*extraordinary},     ## imprecise two-Lorentzian approximation in UV
                {'omega':1.08e15, 'gamma':.15e15, 'sigma':2 + .4*extraordinary}            
                ]
        self.name = "polycrystalline TiO2 (rutile)"
        self.shortname = "TiO2 (rutile)"
        self.where = where
#}}}

class material_SiO2():#{{{
    """ Amorphous SiO2 
    Optical resonances fitted manually to experimental spectra
    Note: Discrete Lorentzian oscillators used predict higher losses in the mid-IR region
    Note2: crystalline SiO2 should be similar, but is slightly birefringent

    From Gunde, M. K.: "Vibrational modes in amorphous silicon dioxide" 
        Physica B: Physics of Condensed Matter, Volume 292, Issue 3-4, p. 286-295.
    """
    def __init__(self, where=None, sigmafactor=1):
        #self.eps = 2.4          ## used if no VIS/UV oscillators defined; may range from 2.3-2.5 [Huber]
        self.eps = 1
        self.pol = [
                {'omega': 447*percm, 'gamma': 49*percm, 'sigma':.923},
                {'omega': 811*percm, 'gamma': 69*percm, 'sigma':.082},
                {'omega':1064*percm, 'gamma': 75*percm, 'sigma':.663},
                {'omega':1165*percm, 'gamma': 80*percm, 'sigma':.058},
                {'omega':1228*percm, 'gamma': 65*percm, 'sigma':.017},
                {'omega': 2.6e15, 'gamma': .5e15, 'sigma': .5},
                {'omega': 2.8e15, 'gamma': .2e15, 'sigma': .55},
                ]
        self.name = "Amorphous silica glass"
        self.shortname = "SiO2 (IR)"
        self.where = where
#}}}

class material_SiC():#{{{
    """ Silicon carbide, or SiC
    THz and optical resonances fitted manually to spectra from databases/experiment
    Note: Discrete Lorentzian oscillators used predict higher losses in the mid-IR region

    http://www.ioffe.rssi.ru/SVA/NSM/Semicond/SiC/basic.html (gives phonon at 25.0 THz?)
    """
    def __init__(self, where=None, sigmafactor=1):
        #self.eps = 2.4          ## used if no VIS/UV oscillators defined; may range from 2.3-2.5 [Huber]
        self.eps = 1
        self.pol = [
                {'omega': 23.75e12, 'gamma': .20e12, 'sigma': 3.4},     # optical phonon
                {'omega': 1.75e15, 'gamma': .4e15, 'sigma': 5.4},
                ]
        self.name = "Silicon carbide"
        self.shortname = "SiC"
        self.where = where
#}}}
class material_Si_NIR():#{{{
    """ Undoped silicon, Si for near IR and visible range 
    Could not approximate the infrared absorption exactly with Lorentzians
    http://www.reading.ac.uk/infrared/library/infraredmaterials/ir-infraredmaterials-si.aspx

    """
    def __init__(self, where=None, sigmafactor=1):
        self.eps = 1          ## flat permittivity for MIR-NIR
        self.pol = [
                {'omega': 1.e15, 'gamma': .15e15, 'sigma': 7.5},
                {'omega': 0.83e15, 'gamma': .05e15, 'sigma': 3.0},
                ]
        self.name = "Silicon for 0.4-1.1 microns"
        self.shortname = "Si (NIR-VIS)"
        self.where = where
#}}}
class material_Si_MIR():#{{{
    """ Undoped silicon, Si for middle infrared range (1.5-5 microns)
    Could not approximate the infrared absorption exactly with Lorentzians

    """
    def __init__(self, where=None, sigmafactor=1):
        self.eps = 11.4          ## flat permittivity for MIR-NIR
        self.pol = [
                {'omega': 1.8e13, 'gamma': .5e13, 'sigma': .3e-3},
                ]
        self.name = "Silicon for 1.5-5 microns"
        self.shortname = "Si (MIR)"
        self.where = where
#}}}
class material_InP():#{{{
    """ Indium Phosphide 
    THz and optical resonances fitted manually to experimental spectra
    Note: Discrete Lorentzian oscillators used predict higher losses in the mid-IR region

    THz optical phonon by
    http://www.ewels.info/science/publications/thesis/html/node52.html
    http://www.ioffe.ru/SVA/NSM/Semicond/InP/basic.html
    Borcherds et al., "Phonon dispersion curves in indium phosphide",

    (The SCOUT database reports permittivity crossing zero around 36.6e12 Hz,
    which seems to be either wrong, or based on some conductivity due to doping.)
    According to Domashevskaya2008 [http://www.sciencedirect.com/science/article/pii/S0921510707005156],
    the longitudinal phonon should be found around 380 cm^-1 = 11.4 THz.


    Note: if we do not know the oscillator strength, but can find the frequencies where 
    permittivity crosses zero and also know the high-frequency permittivity far above the 
    oscillator, then the oscillator strength may be estimated using the Lyddane-Sachs-Teller relation 
        eps_b / eps_s = (omega_L / omega_T)**2
    so
        eps_b - eps_s = [(omega_L / omega_T)**2 - 1] * eps_s
    """
    def __init__(self, where=None, sigmafactor=1):
        self.eps = 1
        self.pol = [
                {'omega': 9.15e12, 'gamma': .5e12, 'sigma': 2.6},     # optical phonon
                {'omega': 1.14e15, 'gamma': .25e15, 'sigma': 6.},
                {'omega': 0.772e15, 'gamma': .20e15, 'sigma': 3.},
                ]
        self.name = "Indium Phosphide"
        self.shortname = "InP"
        self.where = where
#}}}
class material_GaAs():#{{{
    """ Gallium arsenide
    THz and optical resonances fitted manually to experimental spectra
    Note: Discrete Lorentzian oscillators used predict higher losses in the mid-IR region

    Note2: Various sources either state it is transmissive around 10 um (= 33 THz), or absorptive:
        http://www.phoenixinfrared.com/materials/galliumarsenide.html
        http://www.iiviinfrared.com/Optical-Materials/gaas.html
        So the MIR absorption may be caused by some impurities; it could not be properly included 
        in the model presented.

    [http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaAs/basic.html]
        Optical phonon at 35 meV = 8.46 THz  
        permittivity in NIR = 10.89 [ioffe/GaAs]

    [http://www.ioffe.rssi.ru/SVA/NSM/Semicond/GaAs/basic.html]
        TO 8.13, LO 8.79 [Strauch1989] -> strenght calculated
        damping unknown
    """
    def __init__(self, where=None, sigmafactor=1):
        self.eps = 1
        self.pol = [
                {'omega': 8.13e12, 'gamma': 1e12, 'sigma': ((8.79e12/8.13e12)**2-1)*10.89},     # optical phonon damping unknown!
                {'omega': 0.725e15, 'gamma': .06e15, 'sigma': 1.89}, ## manual fit
                {'omega': 0.780e15, 'gamma': .1e15, 'sigma': 2.}, ## manual fit
                {'omega': 1.14e15, 'gamma': .30e15, 'sigma': 6.}, ## manual fit
                ]
        self.name = "Gallium arsenide"
        self.shortname = "GaAs"
        self.where = where
#}}}
class material_Al():#{{{
    """ Drude-Lorentz model for Aluminium

    This model defines aluminium with Drude component and several resonances in 
    near-IR and visible range.
        plasma frequency = 3640 THz [see e. g. Pendry1999] and 
        losses gamma=25e12 (roughly electron scattering frequency)
    """
    def __init__(self, where=None, resistivity=0., eps=0.):
        self.eps = 1. 
        omega0 = 1e-20           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        omega_p = 3640e12       ## plasma frequency of aluminium
        #sigmafactor = 1
        self.pol = [
                #{'omega': 1e6*c*1e-20,   'gamma': 1e6*c*0.037908,   'sigma': 7.6347e+41},
                {'omega': omega0,   'gamma': 1e6*c*0.037908,   'sigma': 7.6347e+41 * 1e-20**2 * (1e6*c)**2 / omega0**2},
                {'omega': 1e6*c*0.13066, 'gamma': 1e6*c*0.26858,    'sigma': 1941},
                {'omega': 1e6*c*1.2453,  'gamma': 1e6*c*0.25165,    'sigma': 4.7065},
                {'omega': 1e6*c*1.4583,  'gamma': 1e6*c*1.0897,     'sigma': 11.396},
                {'omega': 1e6*c*2.8012,  'gamma': 1e6*c*2.7278,     'sigma': 0.55813}]
                #{'omega':1.39e9, 'gamma': 1.0e10, 'sigma':2.6e9}, 
                #{'omega': omega0, 'gamma': 25e12, 'sigma': sigmafactor * omega_p**2 / omega0**2}, # (Lorentz) model for Al
                #{'omega': omega0, 'gamma': 1e12, 'sigma': sigmafactor * omega_p**2 / omega0**2}, # stable? model at THz
                #{'omega':1.39e9, 'gamma': 1e10, 'sigma':1e10}, # Original value, use sigmafactor=1
                #]
        self.name = "Aluminium"
        self.shortname = "Al"
        self.where = where
#}}}

class material_Au(): #{{{
    """ Drude-Lorentz model for gold """
    def __init__(self, where=None, resistivity=0., eps=0.):
        self.eps = 1. 
        omega0 = 1e6*c*1e-20           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        self.pol = [
                {'omega': omega0, 'gamma': 1e6*c*0.042747, 'sigma': 4.0314e+41 * 1e-20**2 * (1e6*c)**2 / omega0**2},
                {'omega':1e6*c*0.33472, 'gamma':1e6*c*0.19438, 'sigma':11.363}, ## sum of Lorentzians = 17.86
                {'omega':1e6*c*0.66944, 'gamma':1e6*c*0.27826, 'sigma':1.1836},
                {'omega':1e6*c*2.3947 , 'gamma':1e6*c*0.7017 , 'sigma':0.65677},
                {'omega':1e6*c*3.4714 , 'gamma':1e6*c*2.0115 , 'sigma':2.6455},
                {'omega':1e6*c*10.743 , 'gamma':1e6*c*1.7857 , 'sigma':2.0148},
                ]
        self.name = "Gold"
        self.shortname = "Au"
        self.where = where
#}}}
class material_Ag():#{{{
    """ Drude-Lorentz model for silver """
    def __init__(self, where=None, resistivity=0., eps=0.):
        self.eps = 1. 
        omega0 = 1e6*c*1e-20           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        self.pol = [
                {'omega': omega0, 'gamma': 1e6*c* 0.038715 , 'sigma': 4.4625e+41 * 1e-20**2 * (1e6*c)**2 / omega0**2},  ## Drude term
                {'omega': 1e6*c*0.6581, 'gamma': 1e6*c*3.1343  , 'sigma':7.9247},
                {'omega': 1e6*c*3.6142, 'gamma': 1e6*c*0.36456 , 'sigma':0.50133},
                {'omega': 1e6*c*6.6017, 'gamma': 1e6*c*0.052426, 'sigma':0.013329},
                {'omega': 1e6*c*7.3259, 'gamma': 1e6*c*0.7388  , 'sigma':0.82655},
                {'omega': 1e6*c*16.365, 'gamma': 1e6*c*1.9511  , 'sigma':1.1133},
                ]
        self.name = "Silver (Drude-Lorentz)"
        self.shortname = "Ag"
        self.where = where
#}}}

class material_Ti():#{{{
    """ Drude-Lorentz model for titanium, edited from Rakic et al., Appl. Opt. 1998 """
    def __init__(self, where=None, resistivity=0., eps=0.):
        self.eps = 1. 
        omega0 = 1e6*c*1e-20           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        self.pol = [
                {'omega': omega0,       'gamma': 1e6*c*0.066137, 'sigma': 5.1166e+40 * 1e-20**2 * (1e6*c)**2 / omega0**2},  ## Drude term
                {'omega': 1e6*c*0.62669,'gamma': 1e6*c*1.8357,   'sigma': 79.136},
                {'omega': 1e6*c*1.2461, 'gamma': 1e6*c*2.0309,   'sigma': 8.7496},
                {'omega': 1e6*c*2.0236, 'gamma': 1e6*c*1.3413,   'sigma': 1.5787},
                {'omega': 1e6*c*1.5671, 'gamma': 1e6*c*1.4211,   'sigma': 0.014077},
                ]
        self.name, self.shortname = "Titanium", "Ti"
        self.where = where
#}}}

## -- Obsoleted or experimental -- 
class material_DrudeMetal_old():#{{{
    """ Defines a generic metal with a Drude model

    In Drude model, the electrons are expected to move freely with some inertia, which
    together with their density determines the plasma frequency omega_p.
    Losses are introduced by letting them to randomly scatter with the frequency gamma.
    The relative permittivity eps_r(omega) of metal can then be calculated as:

        eps_r(omega) = 1 - omega_p**2 / (omega**2 - i*omega*gamma), 
        
    Low-frequency limit of permittivity (eps = eps' + 1j * eps''):
            eps'   --->     - omega_p**2/gamma**2,  
            eps''  --->     omega_p / (omega * gamma),     
    and the high-frequency limit of permittivity:
            eps'   --->     1 - omega**2/omega_p**2         
            eps''  --->     1/omega**3                      
    """

    #def __init__(self, where=None, lfconductivity=15e6, f_c=1e14): XXX
    def __init__(self, where=None, lfconductivity=15e3, f_c=1e14):
        """
        At microwave frequencies (where omega << gamma), the most interesting property of the metal 
        is its conductivity sigma(omega), which may be approximated by a nearly constant real value:
            
            conductivity = epsilon/1j * frequency * eps0 * 2pi

        Therefore, for convenience, we allow the user to specify the low-frequency conductivity sigma0 and compute
        the scattering frequency as:
        """
        #self.gamma = omega_p**2 * epsilon_0 / lfconductivity
        """
        Values similar to those of aluminium are used by default:
            dielectric part of permittivity: 1.0, 
            plasma frequency = 3640 THz
            low-frequency conductivity: 40e6 [S/m]
        """
        ## Design a new metallic model
        ## We need to put the scattering frequency below fc, so that Re(eps) is not constant at f_c
        self.gamma = .5 * f_c           

        ## The virtual plasma frequency is now determined by lfconductivity
        f_p = (self.gamma * lfconductivity / (2*pi) / epsilon_0)**.5       
        self.eps = (f_p/f_c)**2  ## add such an epsilon value, that just shifts the permittivity to be positive at f_c
        #self.eps = 1 ##XXX
        print "F_C", f_c
        print "GAMMA", self.gamma
        print "F_P", f_p
        print "LFC", lfconductivity
        print "eps", self.eps
        print


        ## Feed MEEP with a (fake) Drude model
        f_0 = 1e5           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        self.pol = [
                {'omega': f_0, 'gamma': self.gamma, 'sigma': f_p**2 / f_0**2}, # (Lorentz) model
                ]
        self.name = "Drude metal for <%.2g Hz" % f_c
        self.where = where
#}}}
class material_testsnom():#{{{
    """
    This material is just for testing the s-SNOM response 
    """
    def __init__(self, where=None):
        self.eps = 1.5
        self.pol = [
                #{'omega': 1e10, 'gamma':1e11, 'sigma':100},
                #{'omega':500e12, 'gamma':200e12, 'sigma':1},
                {'omega': 2e13, 'gamma':3e12, 'sigma':42},  
                ]
        self.name = ""
        self.where = where
#}}}
class material_NbN_03K():#{{{  
    """
    Niobium nitride -- low-temperature type-II superconductor (Tk ~ 15 K ?)
    """
    def __init__(self, where=None):
        #self.eps = 100.
        self.eps = 3.
        self.pol = [
                {'omega':2.3e12, 'gamma':.4e12*1.000, 'sigma':5},
                ]
        self.name = "Niobium nitride (T = 3 K)"
        self.where = where
#}}}

class material_AuL():#{{{
    """ Drude-Lorentz model for Gold """
    def __init__(self, where=None, resistivity=0., eps=0.):
        #self.eps = 1. 
        self.eps = 1. 
        omega0 = 1e6*c*1e-20           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        self.pol = [
                {'omega': omega0,	'gamma': 1e6*c*0.042747, 'sigma': 4.0314e+39 * 1e-20**2 * (1e6*c)**2 / omega0**2},
                {'omega':1e6*c*0.33472, 'gamma':1e6*c*0.19438, 'sigma':11.363}, ## sum of Lorentzians = 17.86
                {'omega':1e6*c*0.66944, 'gamma':1e6*c*0.27826, 'sigma':1.1836},
                {'omega':1e6*c*2.3947 , 'gamma':1e6*c*0.7017 , 'sigma': 0.65677},
                {'omega':1e6*c*3.4714 , 'gamma':1e6*c*2.0115 , 'sigma': 2.6455},
                {'omega':1e6*c*10.743 , 'gamma':1e6*c*1.7857 , 'sigma': 2.0148},
                ]
        self.name = "Gold (Drude-Lorentz)"
        self.shortname = "Au"
        self.where = where
#}}}
class material_Au2():#{{{
    """ Drude-Lorentz model for Gold """
    def __init__(self, where=None, resistivity=0., eps=0.):
        #self.eps = 1. 
        self.eps = 18.86 
        omega0 = 1e6*c*1e-20           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        self.pol = [
                {'omega': omega0,	'gamma': 1e6*c*0.042747, 'sigma': 4.0314e+41 * 1e-20**2 * (1e6*c)**2 / omega0**2},
                ]
        self.name = "Gold (Drude-Lorentz)"
        self.shortname = "Au"
        self.where = where
#}}}
class material_Au3():#{{{
    """ Drude-Lorentz model for Gold """
    def __init__(self, where=None, resistivity=0., eps=0.):
        #self.eps = 1. 
        self.eps = 18.86 
        omega0 = 1e6*c*1e-20           ## arbitrary low frequency that makes Lorentz model behave as Drude model
        self.pol = [
                {'omega': omega0,	'gamma': 0, 'sigma': 4.0314e+41 * 1e-20**2 * (1e6*c)**2 / omega0**2},
                ]
        self.name = "Gold (lossy Drude)"
        self.shortname = "Au"
        self.where = where
#}}}

class material_AuL():#{{{
    """ Drude-Lorentz model for Gold - new prototype
    
    
    """
    def __init__(self, where=None, resistivity=0., eps=0.):
        #self.eps = 1. 
        self.eps = 1. 
        self.pol = [
                #{'omega': omega0,	    'gamma': 1e6*c*0.042747, 'sigma': 4.0314e+39 * 1e-20**2 * (1e6*c)**2 / omega0**2},
                {'omega':1e6*c*0.33472, 'gamma':1e6*c*0.19438, 'sigma':11.363}, ## sum of Lorentzians = 17.86
                {'omega':1e6*c*0.66944, 'gamma':1e6*c*0.27826, 'sigma':1.1836},
                {'omega':1e6*c*2.3947 , 'gamma':1e6*c*0.7017 , 'sigma': 0.65677},
                {'omega':1e6*c*3.4714 , 'gamma':1e6*c*2.0115 , 'sigma': 2.6455},
                {'omega':1e6*c*10.743 , 'gamma':1e6*c*1.7857 , 'sigma': 2.0148},
                ]

        #self.Drude_omegap = 4.0314e+39 * 1e-20**2 * (1e6*c)**2)**.5, 
        #self.Drude_gamma = 1e6*c*0.042747

        self.name = "Gold (New Drude-Lorentz)"
        self.shortname = "Au"
        self.where = where
#}}}

class material_SiO2_stability_test():#{{{
    """ Amorphous SiO2  for mid_infrared
    Optical resonances fitted manually to experimental spectra
    Note: Discrete Lorentzian oscillators used predict higher losses in the mid-IR region
    Note2: crystalline SiO2 should be similar, but is slightly birefringent

    From Gunde, M. K.: "Vibrational modes in amorphous silicon dioxide" 
        Physica B: Physics of Condensed Matter, Volume 292, Issue 3-4, p. 286-295.
    """
    def __init__(self, where=None, sigmafactor=1):
        #self.eps = 2.4          ## used if no VIS/UV oscillators defined; may range from 2.3-2.5 [Huber]
        #self.eps = 1+.5+.55+ .082+ .663+ .058+ .017

        ## S = 0.5783  eps(f_c)=1                   UNSTAB, ok
        ## S = 0.5763  eps(f_c)=1                   stab      
        ## S = 0.5763  eps(f_c)=.99                 stab, WEIRD!
        ## S = 0.5783  eps(f_c)=.95                 UNSTAB, ok
        ## S = 0.5763  eps(f_c)=.9989+0.015j        UNSTAB, WEIRD!
        #M S = 0.5763  eps(f_c)=.9989+0.015j        UNSTAB, WEIRD!

        omega0 = 3.0e-10
        omega = 3.0e13
        self.eps = 1.42

        self.pol = [
                #{'omega': 447*percm, 'gamma': 49*percm, 'sigma':.923},
                {'omega': omega0, 'gamma': 9e12, 'sigma':.523 * (omega/omega0)**2},
                ]
        self.name = "Amorphous silica glass (SiO2) for IR range"
        self.shortname = "SiO2 (IR)"
        self.where = where
#}}}

class material():#{{{
    """ Base class for materials """
    def __init__(self, where=None):
        self.eps = 1
        self.pol = [
                ]
        self.name = "default material"
        self.shortname = "default"
        self.where = where

#}}}
