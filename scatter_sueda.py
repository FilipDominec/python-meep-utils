#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" An example python-meep simulation of a dielectric sphere scattering a broadband impulse, 
illustrating the use of the convenient functions provided by meep_utils.py 
(c) 2014 Filip Dominec, see http://f.dominec.eu/meep for more information """

import time, sys, os
import numpy as np
from scipy.constants import c, epsilon_0, mu_0

import meep_utils, meep_materials, metamaterial_models
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep

class MetalInterface(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=15e-15, resolution=4e-9, cells=1, cell_size=50e-9, padding=20e-9):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation

        ## Constant parameters for the simulation
        self.simulation_name = "metalinterface"    
        self.src_freq, self.src_width = 1000e12, 4000e12  # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (100e12, 2000e12)    # Which frequencies will be saved to disk
        self.pml_thickness = 20e-9

        self.size_x = resolution*2 
        self.size_y = resolution*2
        self.size_z = cells*cell_size + 4*padding + 2*self.pml_thickness
        self.monitor_z1, self.monitor_z2 = (-(cell_size*cells/2)-padding, (cell_size*cells/2)+padding)
        self.cellcenters = np.arange((1-cells)*cell_size/2, cells*cell_size/2, cell_size)

        self.register_locals(locals())          ## Remember the parameters

        ## Define materials
        self.materials = []  
        if 'Au' in comment:           
             self.materials += [meep_materials.material_Au(where=self.where_metal),
                    meep_materials.material_Ag(where=None)]
        else:
             self.materials += [meep_materials.material_Au(where=None), # meep_materials.material_Au(where=None),
                    meep_materials.material_dielectric(eps=10, where=self.where_metal)]

        for m in self.materials: 
            self.fix_material_stability(m, f_c=5e15, verbose=0) ## rm all osc above the first one, to optimize for speed 

        ## Test the validity of the model
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_metal(self, r):
        if in_zslab(r, d=self.cell_size/2, cz=0):
            return self.return_value             # (do not change this line)
        return 0
#}}}

# Model selection
sim_param, model_param = meep_utils.process_param(sys.argv[1:])
#model = metamaterial_models.SphereArray(**model_param)
model = MetalInterface(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume, structure and the fields according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes=meep.Z)
f = meep.fields(s)
f.use_bloch(meep.X, sim_param.get('Kx', 0) / (-2*np.pi)) # (any transversal component of k-vector is allowed)
f.use_bloch(meep.Y, sim_param.get('Ky',.0) / (-2*np.pi))

# Add the field source (see meep_utils for an example of how an arbitrary source waveform is defined)
if not sim_param['frequency_domain']:       ## (temporal source shape)
    #src_time_type = meep.band_src_time(model.src_freq/c, model.src_width/c, model.simtime*c/1.1)
    src_time_type = meep.gaussian_src_time(model.src_freq/c, model.src_width/c)
else:
    src_time_type = meep.continuous_src_time(sim_param['frequency']/c)
srcvolume = meep.volume(                    ## (spatial source shape)
        meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness),
        meep.vec( model.size_x/2,  model.size_y/2, -model.size_z/2+model.pml_thickness))
f.add_volume_source(meep.Ex, src_time_type, srcvolume)

## Define monitors planes and visualisation output
monitor_options = {'size_x':model.size_x, 'size_y':model.size_y, 'Kx':sim_param.get('Kx', 0), 'Ky':sim_param.get('Ky', 0)}
monitor1_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z1, **monitor_options)
monitor1_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z1, **monitor_options)
monitor2_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z2, **monitor_options)
monitor2_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z2, **monitor_options)

slices = []
slices += [meep_utils.Slice(model=model, field=f, components=(meep.Dielectric), at_t=0, name='EPS')]
#slices += [meep_utils.Slice(model=model, field=f, components=(meep.Ex), at_x=0, name='FieldEvolution', min_timestep=1e-12)]
slices += [meep_utils.Slice(model=model, field=f, components=(meep.Ex, meep.Ey, meep.Ez), at_t=np.inf, name='SnapshotE')]
slices += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_x=0, at_t=np.inf, 
    name=('At%.3eHz'%sim_param['frequency']) if sim_param['frequency_domain'] else '', outputpng=True, outputvtk=False)]

## Run the FDTD simulation or the frequency-domain solver
if not sim_param['frequency_domain']:       ## time-domain computation
    f.step(); timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    while (f.time()/c < model.simtime):     # timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        #print f.get_field(meep.Ex, meep.vec(0,0,0))
        for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
        for slice_ in slices: slice_.poll(f.time()/c)
    for slice_ in slices: slice_.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    f.step()
    f.solve_cw(sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab']) 
    for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
    for slice_ in slices: slice_.finalize()
    meep_utils.notify(model.simulation_name)

## Get the reflection and transmission of the structure
if meep.my_rank() == 0:
    freq, s11, s12, headerstring = meep_utils.get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, 
            frequency_domain=sim_param['frequency_domain'], frequency=sim_param['frequency'], 
            intf=getattr(model, 'interesting_frequencies', [0, model.src_freq+model.src_width]),
            pad_zeros=1.0, Kx=sim_param.get('Ky', 0), Ky=sim_param.get('Ky', 0))

    meep_utils.savetxt(fname=model.simulation_name+".dat", fmt="%.6e",
            X=zip(freq, np.abs(s11), np.angle(s11), np.abs(s12), np.angle(s12)), 
            header=model.parameterstring + meep_utils.sim_param_string(sim_param) + headerstring)

    with open("./last_simulation_name.dat", "w") as outfile: outfile.write(model.simulation_name) 

meep.all_wait()         # Wait until all file operations are finished
