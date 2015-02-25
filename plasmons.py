#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" run_sim.py - an example python-meep simulation of a thin metal sheet with two apertures, coupling
the back-incident light to surface plasmons which then show interference pattern.
Illustrates the use of the convenient functions provided by meep_utils.py 
(c) 2015 Filip Dominec, see http://fzu.cz/~dominecf/meep/ for more information """

import numpy as np
import time, sys, os
import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep
c = 2.997e8

sim_param, model_param = meep_utils.process_param(sys.argv[1:])
class plasmonyp_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-15, resolution=100e-9, size_x=16e-6, size_y=5e-6, size_z=5e-6,
            metalthick=.5e-6, apdisty=0e-6, apdistx=14e-6, aprad=.1e-6, monzd=123.456e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "PlasmonsYPb"    
        monzd=size_z

        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        substrate_z = size_x / 3
        self.pml_thickness = 1e-6
        self.monitor_z1, self.monitor_z2 = (-(monzd)/2, (monzd)/2)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 400e12, 800e12     # [Hz] (note: gaussian source ends at t=10/srcWidth)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk
        self.Kx = 0; self.Ky = 0; self.padding=0
        self.size_x = size_x 
        self.size_y = size_y
        self.size_z = size_z

        ## Define materials
        f_c = c / np.pi/self.resolution/meep_utils.meep.use_Courant()

        self.materials   = [meep_materials.material_Au(where=self.where_metal)]  
        #self.materials   = [meep_materials.material_DrudeMetal(lfconductivity=1e8, f_c=.2*f_c, where = self.where_metal)]  
        self.materials   += [meep_materials.material_dielectric(where=self.where_diel, eps=2.)]  
        self.TestMaterials()
        self.materials   += [meep_materials.material_Au(where=None)]  
        meep_utils.plot_eps(self.materials, mark_freq=[f_c])

    # each material has one callback, used by all its polarizabilities (thus materials should never overlap)
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

#}}}

# Model selection
model = plasmonyp_model(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume and structure according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
#s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes=meep.Z)
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes="All")

## Create fields with Bloch-periodic boundaries (any transversal component of k-vector is allowed, but may not radiate)
f = meep.fields(s)

## Add a source of the plane wave (see meep_utils for definition of arbitrary source shape)
if not sim_param['frequency_domain']:           ## Select the source dependence on time
    #src_time_type = meep.band_src_time(model.srcFreq/c, model.srcWidth/c, model.simtime*c/1.1)
    #src_time_type = meep.gaussian_src_time(model.srcFreq/c, model.srcWidth/c)
    src_time_type = meep.continuous_src_time(model.srcFreq/c)
else:
    src_time_type = meep.continuous_src_time(sim_param['frequency']/c)
srcvolume = meep.volume( 
        meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness),
        meep.vec(+model.size_x/2, +model.size_y/2, -model.size_z/2+model.pml_thickness))
f.add_volume_source(meep.Ex, src_time_type, srcvolume)


## Define visualisation output
slice_makers =  [meep_utils.Slice(model=model, field=f, components=(meep.Dielectric), at_t=0, name='EPS')]
slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ez, at_y=0, min_timestep=.3e-15, outputgif=True, name='ParallelCut')]
slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ez, at_z=model.metalthick/2+model.resolution, min_timestep=.3e-15, outputgif=True, name='PerpendicularCut')]
slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_t=100e-15)]

if not sim_param['frequency_domain']:       ## time-domain computation
    f.step()
    dt = (f.time()/c)
    meep_utils.lorentzian_unstable_check_new(model, dt)
    timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    while (f.time()/c < model.simtime):                               # timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        meep.master_printf("%e" % abs(f.get_field(meep.Ex, meep.vec(model.size_x/4, model.size_y/4, model.size_z/4))))
        for slice_maker in slice_makers: slice_maker.poll(f.time()/c)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    f.step()
    f.solve_cw(sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab']) 
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name)

with open("./last_simulation_name.txt", "w") as outfile: outfile.write(model.simulation_name) 
meep.all_wait()         # Wait until all file operations are finished
