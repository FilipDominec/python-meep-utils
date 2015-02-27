#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" (c) 2014 Filip Dominec, see http://fzu.cz/~dominecf/meep/ for more information """


import numpy as np
import time, sys, os
import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep
c = 2.997e8

sim_param, model_param = meep_utils.process_param(sys.argv[1:])
class HollowCyl_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=30e-9, resolution=5e-3, cells=1, padding=9e-3, radius=33.774e-3, height=122.36e-3 ,  Kx=0, Ky=0): ## XXXheight=122.3642686e-3
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "HollowCyl"    
        
        self.register_locals(locals())          ## Remember the parameters

        ## Obligatory parameters (used in the simulation)
        self.pml_thickness = padding/2
        self.simtime = simtime      # [s]HollowCyl_simtime=3.000e-08_height=3.000e-02
        self.src_freq, self.src_width = 3e9, 10e9     # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (.1e9, 8e9)     # Which frequencies will be saved to disk
        self.size_x = 2*radius+padding*2
        self.size_y = 2*radius+padding*2
        self.size_z = height+padding*2

        ## Define materials
        f_c = c / np.pi/self.resolution/meep_utils.meep.use_Courant() 
        self.materials = []  
        self.materials += [meep_materials.material_DrudeMetal(lfconductivity=1e4, f_c=f_c, gamma_factor=.5, epsplus=0, where=self.where_metal)]  
        #self.materials += [meep_materials.material_dielectric(eps=.5**.5, where = None)]  
        #self.materials += [meep_materials.material_SiO2_stability_test(where = self.where_TiO2)]  

        #meep.use_Courant(3**-.5 - 0.001); print "COURANT =", meep.use_Courant()
        print "model thinks f_c = ", f_c

        meep_utils.plot_eps(self.materials, mark_freq={f_c:'$f_c$'}, plot_conductivity=True) #10e12:' $\\frac{\\gamma}{2\\pi}$ ' , .7e15:'$\\frac{\\omega_p}{2\\pi}$ '
        #eps_fc = meep_utils.analytic_eps(self.materials[0], f_c)
        #lim_tmp = meep.use_Courant()**2 * 3
        #print "analytic eps(f_c) = ", eps_fc, (" [the limit is %f]" % lim_tmp), "******* UNSTABLE ******* " if (eps_fc.real< lim_tmp) else ""

        self.test_materials()

    def where_metal(self, r):
        if not (in_zcyl(r, cx=0, cy=0, rad=self.radius) and in_zslab(r, cz=0, d=self.height)):  # 
        #if  in_zslab(r, cz=0, d=self.radius):
            return self.return_value             # (do not change this line)
        return 0
#}}}


# Model selection
model = HollowCyl_model(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume and structure according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes="All") ## XXX   meep.XY


## Create the fields object, and define the Bloch-periodic boundaries (any transversal component of k-vector is allowed)
field = meep.fields(s)
#field.use_bloch(meep.Z, 1)          ## periodic along the cylinder XXX
# (removing cylinder caps -> making an infinite waveguide with periodic boundary condition, change pml_axes=meep.XY)

# Add the field source (see meep_utils for an example of how an arbitrary source waveform is defined)
if not sim_param['frequency_domain']:           ## Select the source dependence on time
    #src_time_type = meep.band_src_time(-model.src_freq/c, model.src_width/c, model.simtime*c/1.1)
    src_time_type = meep.gaussian_src_time(-model.src_freq/c, model.src_width/c)  ## negative frequency supplied -> e^(+i omega t) convention
else:
    src_time_type = meep.continuous_src_time(-sim_param['frequency']/c) ## TODO check in freq domain that negative frequency is OK, and is it needed?
srcvolume = meep.volume( 
        meep.vec(model.radius*.15, -model.radius*.25, -model.height*.15),
        meep.vec(model.radius*.15, -model.radius*.25, -model.height*.15))
field.add_volume_source(meep.Ez, src_time_type, srcvolume) ## source of oblique polarization - excites both TE and TM modes
field.add_volume_source(meep.Ex, src_time_type, srcvolume)

#monitor_options = {'size_x':model.size_x/4, 'size_y':model.size_y/4, 'Kx':model.Kx, 'Ky':model.Ky}
#monitor1_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.size_z/8, **monitor_options)

slice_makers =  []
#slice_makers += [meep_utils.Slice(model=model, field=field, components=(meep.Ex, meep.Ey, meep.Ez), at_t=3, name="ElectricAtEnd")]
#slice_makers += [meep_utils.Slice(model=model, field=field, components=(meep.Hx, meep.Hy, meep.Hz), at_t=3, name="MagneticAtEnd")]
#slice_makers =  [meep_utils.Slice(model=model, field=field, components=(meep.Dielectric), at_t=0, name='EPS')]

if not sim_param['frequency_domain']:       ## time-domain computation
    field.step()
    dt = (field.time()/c)
    meep_utils.lorentzian_unstable_check_new(model, dt, quit_on_warning=False)
    timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    monitor_point = meep.vec(-model.radius*.5, model.radius*.3, model.height*.3)
    x,y = [], []
    while (field.time()/c < model.simtime):                               # timestepping cycle
        field.step()
        timer.print_progress(field.time()/c)
        if field.time()/c > 30/model.src_width:
            x.append(field.time()/c); 
            y.append(field.get_field(meep.Ex, monitor_point)+field.get_field(meep.Ey, monitor_point)+field.get_field(meep.Ez, monitor_point))
        for slice_maker in slice_makers: slice_maker.poll(field.time()/c)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    field.step()
    field.solve_cw(sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab']) 
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name)

# Get the reflection and transmission of the structure
if meep.my_rank() == 0 and  not sim_param['frequency_domain']:
    x, y = np.array(x), np.array(y)
    ## Plot time-domain data
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10,10))
    plt.plot(x, np.abs(y), color="k", label=u"$|y|$", ls='-')       
    plt.xlabel(u"time [s]"); plt.ylabel(u"amplitude"); plt.yscale('log'); plt.grid()
    plt.legend(prop={'size':10}, loc='upper right').draw_frame(False)
    plt.savefig("td.png", bbox_inches='tight')

    ## 1D FFT with cropping for useful frequencies
    plt.figure(figsize=(15,10))
    freq    = np.fft.fftfreq(len(x), d=(x[1]-x[0]))         # calculate the frequency axis with proper spacing
    yf      = np.fft.fft(y, axis=0) / len(x) * 2*np.pi      # calculate the FFT values
    freq    = np.fft.fftshift(freq)                         # ensures the frequency axis is a growing function
    yf      = np.fft.fftshift(yf) / np.exp(1j*2*np.pi*freq * x[0])   # dtto, and corrects the phase for the case when x[0] != 0
    truncated = np.logical_and(freq>0, freq<model.interesting_frequencies[1])         # (optional) get the frequency range
    (yf, freq) = map(lambda array: array[truncated], (yf, freq))    # (optional) truncate the data points
    plt.plot(freq, np.abs(yf), color="#FF8800", label=u"$y'$", ls='-')                  # (optional) plot amplitude


    ## Convert to polar notation and save
    meep_utils.savetxt(fname=model.simulation_name+".dat", X=zip(freq, abs(yf), meep_utils.get_phase(yf)), fmt="%.8e",
            header=model.parameterstring+"#x-column _frequency [Hz]\n#column ampli\n#column phase\n")


    ## Harminv
    ## TODO switch to harminv_wrapper.py instead
    hi = meep_utils.harminv(x, y, amplitude_prescaling=1e6)
    oscillator_count = len(hi['frequency'])
    #if oscillator_count > 0:
        #plt.scatter(np.abs(hi['frequency']), hi['amplitude'], c=hi['phase'], s=np.abs(hi['quality'])/50 + 2, cmap=plt.cm.hsv, alpha=.2)

    def lorentz(x, f, d, q, A, p, err):
        return A*abs(q)/(1 + abs(q*np.pi*2)*(x-abs(f))**2) / (np.pi*2) * 1e16 # XXX

    freq_fine = np.linspace(0, np.max(freq)*1, 1000)
    sumosc = np.zeros_like(freq_fine)
    print "Harminv frequencies", np.abs(hi['frequency'])
    for osc in range(oscillator_count):
        osc_y = lorentz(freq_fine,   hi['frequency'][osc], hi['decay'][osc], hi['quality'][osc], hi['amplitude'][osc], hi['phase'][osc], hi['error'][osc])
        plt.plot(freq_fine, osc_y, color="#0088FF", label=u"", ls='-', alpha=.3)      # (optional) plot amplitude
        sumosc += osc_y 
    plt.plot(freq_fine, sumosc, color="#0088FF", label=u"$\\Sigma$ osc", ls='-')      # (optional) plot amplitude

    analytic_modes = {}
    from scipy.special import jnyn_zeros
    # For a long-enough cavity, the lines group as such: [Pozar: microwave engineering],  B'01 = 3.832
    # TE111-TE112      TM010-TM011-TM012     TE211-TE212     TM110-TM111+TE011-TM112
    # so TE for p=0 (TExx0) is not allowed
    #    TM for m=0 (TMx0x) is not allowed
    # In the plot of mine, they grou as such:
    # TE101-TE102-TE103    TM000-TM002-TM003        TE201-TE202-TE203       TM101-
    
    freq_correction = (1. - 150e6/3.2e9)
    for p in [0,1,2]:
        for n in range(4):
            S = " "*p
            for m,B in enumerate(jnyn_zeros(n, 5)[0]):
                analytic_modes[freq_correction * c/(2*np.pi) * np.sqrt((B/(model.radius))**2 + (p*np.pi/model.height)**2) ] = ("%s$TM_{%d%d%d}$%s" % (S,n,m+1,p,S))
            for m,B in enumerate(jnyn_zeros(n, 5)[1]):
                if p>0: ## TExx0 modes can not exist [Pozar]
                    analytic_modes[freq_correction * c/(2*np.pi) * np.sqrt((B/(model.radius))**2 + (p*np.pi/model.height)**2) ] = ("   %s$TE_{%d%d%d}$%s   " % (S,n,m+1,p,S))
    print analytic_modes
    meep_utils.annotate_frequency_axis(analytic_modes, label_position_y=1, arrow_length=10, log_y=True)

    ## Finish the plot + save 
    plt.xlim(model.interesting_frequencies); 
    plt.xlabel(u"frequency [Hz]"); 
    plt.ylabel(u"amplitude excited by a pulse"); 
    plt.ylim((1e-3, 1e4))
    plt.yscale('log')
    plt.grid()
    plt.legend(prop={'size':10}, loc='upper right').draw_frame(False)
    plt.savefig("%s.png" % model.simulation_name, bbox_inches='tight')


    #meep_utils.savetxt(freq=freq, s11=s11, s12=s12, model=model)
    #with open("./last_simulation_name.txt", "w") as outfile: outfile.write(model.simulation_name) 
    #import effparam        # process effective parameters for metamaterials

meep.all_wait()         # Wait until all file operations are finished
