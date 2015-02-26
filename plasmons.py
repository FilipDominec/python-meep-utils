#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" run_sim.py - an example python-meep simulation of a thin metal sheet with an aperture, coupling
the back-incident light to surface plasmons. If the film is surrounded by two media with similar index
of refraction, circular interference pattern can be observed between the symmetric and antisymmetric 
plasmon modes. 
A different (hyperbolic) interference pattern can be obtained when the plasmons are radiated by two
holes, even  if there is a mismatch between the refractive indices of the substrate and superstrate.
(c) 2015 Filip Dominec, see http://f.dominec.eu/meep for more information """

import time, sys, os
import numpy as np
from scipy.constants import c, epsilon_0, mu_0

import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep

class PlasmonFilm_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-15, resolution=100e-9, size_x=16e-6, size_y=5e-6, size_z=5e-6,
            metalthick=.5e-6, apdisty=0e-6, apdistx=14e-6, aprad=.1e-6, monzd=123.456e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation

        ## Constant parameters for the simulation
        self.simulation_name = "PlasmonsFilm"    
        self.src_freq, self.src_width = 400e12, 800e12     # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk
        self.pml_thickness = 1e-6

        self.size_x = size_x 
        self.size_y = size_y
        self.size_z = size_z
        substrate_z = size_x / 3
        self.simtime = simtime      # [s]
        self.Kx = 0; self.Ky = 0; self.padding=0
        self.register_locals(locals())          ## Remember the parameters
        ## Define materials
        f_c = c / np.pi/self.resolution/meep_utils.meep.use_Courant()

        self.materials   = [meep_materials.material_Au(where=self.where_metal)]  
        #self.materials   = [meep_materials.material_DrudeMetal(lfconductivity=1e8, f_c=.2*f_c, where = self.where_metal)]  
        self.materials   += [meep_materials.material_dielectric(where=self.where_diel, eps=2.)]  
        #self.TestMaterials()
        #self.materials   += [meep_materials.material_Au(where=None)]  
        meep_utils.plot_eps(self.materials, mark_freq=[f_c])

        ## Test the validity of the model
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_metal(self, r):
        if in_zslab(r, cz=0, d=self.metalthick) and not \
                (in_zcyl(r, cx=-self.apdistx/2,  cy=-self.apdisty/2, rad=self.aprad) or \
                in_zcyl(r,  cx=1000+ self.apdistx/2,  cy= self.apdisty/2, rad=self.aprad)): ## XXX second aperture disabled!
            return self.return_value             # (do not change this line)
        return 0
    def where_diel(self, r):
        if r.z() < 0:
            return self.return_value             # (do not change this line)
        return 0
 # TODO |r|   r phase
#}}}

# Model selection
sim_param, model_param = meep_utils.process_param(sys.argv[1:])
model = PlasmonFilm_model(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume, structure and the fields according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes="All")

## Create fields with Bloch-periodic boundaries 
f = meep.fields(s)

# Add the field source (see meep_utils for an example of how an arbitrary source waveform is defined)
src_time_type = meep.continuous_src_time(model.src_freq/c)
srcvolume = meep.volume( 
        meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness),
        meep.vec(+model.size_x/2, +model.size_y/2, -model.size_z/2+model.pml_thickness))
f.add_volume_source(meep.Ex, src_time_type, srcvolume)


## Define visualisation output
slices =  [meep_utils.Slice(model=model, field=f, components=(meep.Dielectric), at_t=0, name='EPS')]
slices += [meep_utils.Slice(model=model, field=f, components=meep.Ez, at_y=0, min_timestep=.3e-15, outputgif=True, name='ParallelCut')]
slices += [meep_utils.Slice(model=model, field=f, components=meep.Ez, at_z=model.metalthick/2+model.resolution, min_timestep=.3e-15, outputgif=True, name='PerpendicularCut')]
slices += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_t=100e-15)]

if not sim_param['frequency_domain']:       ## time-domain computation
    f.step(); timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    while (f.time()/c < model.simtime):     # timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        #meep.master_printf("%e" % abs(f.get_field(meep.Ex, meep.vec(model.size_x/4, model.size_y/4, model.size_z/4))))
        for slice_ in slices: slice_.poll(f.time()/c)
    for slice_ in slices: slice_.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    f.step()
    f.solve_cw(sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab']) 
    for slice_maker in slices: slice_maker.finalize()
    meep_utils.notify(model.simulation_name)


if meep.my_rank() == 0:
    with open("./last_simulation_name.dat", "w") as outfile: outfile.write(model.simulation_name) 
meep.all_wait()         # Wait until all file operations are finished
