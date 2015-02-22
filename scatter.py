#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" An example python-meep simulation of a dielectric sphere scattering a broadband impulse, 
illustrating the use of the convenient functions provided by meep_utils.py 
(c) 2014 Filip Dominec, see http://f.dominec.eu/meep for more information """


import numpy as np
import time, sys, os
import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep
c = 2.997e8

sim_param, model_param = meep_utils.process_param(sys.argv[1:])
class SphereArray_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=2e-12, resolution=2e-6, cells=1, cell_size=50e-6, padding=20e-6, Kx=0, Ky=0, 
            radius=13e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SphereArray"    
        
        ## Constants for the simulation
        self.pml_thickness = 20e-6
        #self.monitor_z1, self.monitor_z2 = (-(cell_size*cells/2)-padding, (cell_size*cells/2)+padding)
        self.monitor_z1, self.monitor_z2 = ((-cells*cell_size/2)-padding, (cells*cell_size/2)+padding)
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 1000e9, 2000e9     # [Hz] (note: gaussian source ends at t=10/srcWidth)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk

        self.layer_thickness = (self.monitor_z2 - self.monitor_z1 - 2*padding) # TODO layer_thickness identical with cell_size
        #if "interesting_frequencies" in dir(model): interest_freq = self.interesting_frequencies 
        #else: 
        #interest_freq = (0, self.srcFreq+self.srcWidth) # TODO remove
        #self.plot_freq_min =  interest_freq[0]
        #self.plot_freq_max =  interest_freq[1]
        self.simulation_orig_name = self.simulation_name # TODO remove

        self.register_locals(locals())          ## Remember the parameters

        self.size_x = self.cell_size 
        self.size_y = self.cell_size
        self.size_z = self.cells*self.cell_size + 2*self.pml_thickness + 4*self.padding

        ## Define materials
        self.materials = []  
        #self.materials += [meep_materials.material_Au(where=self.where_sphere)]  
        self.materials += [meep_materials.material_dielectric(where=self.where_sphere, eps=100)]  
        #for x in [.01, .1, 1, 10]: ## TODO remove
            #self.materials += [make_stable_Drude(meep_materials.material_Au(where = None), f_c=f_c*x)]  
            #self.materials += [make_stable_Drude(meep_materials.material_AuL(where = None), f_c=f_c*x)]  
            #self.materials += [meep_materials.material_DrudeMetal(where = self.where_TiO2, lfconductivity=2.5e6*2*np.pi, f_c=f_c)]  

        for n, material in enumerate(self.materials): self.fix_material_stability(material)
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_block=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_sphere(self, r):
        if  in_sphere(r, cx=0, cy=0, cz=0, rad=self.radius):
        #if  in_zslab(r, cz=0, d=self.radius):
            return self.return_value             # (do not change this line)
        return 0
        #print "analytic eps(f_c) = ", eps_fc, (" [the limit is %f]" % lim_tmp), "******* UNSTABLE ******* " if (eps_fc.real< lim_tmp) else ""
        #eps_fc = meep_utils.analytic_eps(self.materials[0], f_c)
        #lim_tmp = meep.use_Courant()**2 * 3
#}}}



# Model selection
model = SphereArray_model(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume and structure according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes=meep.Z)


## Create the fields object
f = meep.fields(s)
# Define the Bloch-periodic boundaries (any transversal component of k-vector is allowed)
f.use_bloch(meep.X, -model.Kx/(2*np.pi))
f.use_bloch(meep.Y, -model.Ky/(2*np.pi))

# Add the field source (see meep_utils for an example of how an arbitrary source waveform is defined)
if not sim_param['frequency_domain']:           ## Select the source dependence on time
    #src_time_type = meep.band_src_time(model.srcFreq/c, model.srcWidth/c, model.simtime*c/1.1)
    src_time_type = meep.gaussian_src_time(model.srcFreq/c, model.srcWidth/c)
else:
    src_time_type = meep.continuous_src_time(sim_param['frequency']/c)
srcvolume = meep.volume( 
        meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness),
        meep.vec( model.size_x/2,  model.size_y/2, -model.size_z/2+model.pml_thickness))
f.add_volume_source(meep.Ex, src_time_type, srcvolume)


## Define monitors planes and visualisation output
monitor_options = {'size_x':model.size_x, 'size_y':model.size_y, 'Kx':model.Kx, 'Ky':model.Ky}
monitor1_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z1, **monitor_options)
monitor1_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z1, **monitor_options)
monitor2_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z2, **monitor_options)
monitor2_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z2, **monitor_options)

slice_makers =  [meep_utils.Slice(model=model, field=f, components=(meep.Dielectric), at_t=0, name='EPS')]
slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_x=0, min_timestep=.05e-12)]

if not sim_param['frequency_domain']:       ## time-domain computation
    f.step()
    timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    while (f.time()/c < model.simtime):                               # timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
        for slice_maker in slice_makers: slice_maker.poll(f.time()/c)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    f.step()
    f.solve_cw(sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab']) 
    for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name)

## Get the reflection and transmission of the structure
if meep.my_rank() == 0:

    freq, s11, s12 = meep_utils.get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, 
            frequency_domain=sim_param['frequency_domain'], frequency=sim_param['frequency'], 
            maxf=model.srcFreq+model.srcWidth, pad_zeros=1.0, Kx=model.Kx, Ky=model.Ky)

    meep_utils.savetxt(fname=model.simulation_name+".dat", X=zip(freq, np.abs(s11), np.angle(s11), np.abs(s12), np.angle(s12)), header=model.parameterstring)

    with open("./last_simulation_name.dat", "w") as outfile: outfile.write(model.simulation_name) 
    import effparam        # process effective parameters for metamaterials

meep.all_wait()         # Wait until all file operations are finished
