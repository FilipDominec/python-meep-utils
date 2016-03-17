#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" (c) 2014 Filip Dominec, see http://fzu.cz/~dominecf/meep/ for more information """

import time, sys, os
import numpy as np
from scipy.constants import c, epsilon_0, mu_0

import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep

model_param = meep_utils.process_param(sys.argv[1:])
class HollowCyl_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=30e-9, resolution=5e-3, cellnumber=1, padding=9e-3, radius=33.774e-3, height=122.36e-3 ,  Kx=0, Ky=0, **other_args):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "HollowCyl"    
        
        self.register_locals(locals(), other_args)          ## Remember the parameters

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
        au = meep_materials.material_Au(where=self.where_metal)
        self.fix_material_stability(au, verbose=0)
        self.materials.append(au)
        #self.materials += [meep_materials.material_DrudeMetal(lfconductivity=1e4, f_c=f_c, gamma_factor=.5, epsplus=0, where=self.where_metal)]  

        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    def where_metal(self, r):
        if not (in_zcyl(r, cx=0, cy=0, rad=self.radius) and in_zslab(r, cz=0, d=self.height)):  # 
        #if  in_zslab(r, cz=0, d=self.radius):
            return self.return_value             # (do not change this line)
        return 0
#}}}

# Model selection
model_param = meep_utils.process_param(sys.argv[1:])
model = HollowCyl_model(**model_param)

## Initialize volume, structure and the fields according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, pml_axes="All") ## XXX   meep.XY
f = meep.fields(s)
# Define the Bloch-periodic boundaries (any transversal component of k-vector is allowed)
#f.use_bloch(meep.Z, 1)          ## periodic along the cylinder XXX
# (removing cylinder caps -> making an infinite waveguide with periodic boundary condition, change pml_axes=meep.XY)

# Add the field source (see meep_utils for an example of how an arbitrary source waveform is defined)
if not getattr(model, 'frequency', None):           ## Select the source dependence on time
    #src_time_type = meep.band_src_time(-model.src_freq/c, model.src_width/c, model.simtime*c/1.1)
    src_time_type = meep.gaussian_src_time(-model.src_freq/c, model.src_width/c)  ## negative frequency supplied -> e^(+i omega t) convention
else:
    src_time_type = meep.continuous_src_time(-getattr(model, 'frequency', None)/c) ## TODO check in freq domain that negative frequency is OK, and is it needed?
srcvolume = meep.volume( 
        meep.vec(model.radius*.15, -model.radius*.25, -model.height*.15),
        meep.vec(model.radius*.15, -model.radius*.25, -model.height*.15))
f.add_volume_source(meep.Ez, src_time_type, srcvolume) ## source of oblique polarization - excites both TE and TM modes
f.add_volume_source(meep.Ex, src_time_type, srcvolume)

slices = []
if not "noepssnapshot" in model.comment:
    slices += [meep_utils.Slice(model=model, field=f, components=(meep.Dielectric), at_t=0, name='EPS')]
if "narrowfreq-snapshots" in model.comment:
    slices += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_y=0, at_t=np.inf,
            name=('At%.3eHz'%getattr(model, 'frequency', None)) if getattr(model, 'frequency', None) else '',
            outputpng=True, outputvtk=False)]
if "fieldevolution" in model.comment: 
    slices += [meep_utils.Slice(model=model, field=f, components=(meep.Ex), at_x=0, name='FieldEvolution', 
        min_timestep=.1/model.src_freq, outputgif=True, outputvtk=True)]
if "snapshote" in model.comment:
    slices += [meep_utils.Slice(model=model, field=f, components=(meep.Ex, meep.Ey, meep.Ez), at_t=np.inf, name='SnapshotE')]


if not getattr(model, 'frequency', None):       ## time-domain computation
    f.step(); timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    monitor_point = meep.vec(-model.radius*.5, model.radius*.3, model.height*.3)
    x,y = [], []
    while (f.time()/c < model.simtime):     # timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        if f.time()/c > 30/model.src_width:
            x.append(f.time()/c); 
            y.append(f.get_field(meep.Ex, monitor_point)+f.get_field(meep.Ey, monitor_point)+f.get_field(meep.Ez, monitor_point))
        for slice_ in slices: slice_.poll(f.time()/c)
    for slice_ in slices: slice_.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    f.solve_cw(getattr(model, 'MaxTol',0.001), getattr(model, 'MaxIter', 5000), getattr(model, 'BiCGStab', 8)) 
    for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
    for slice_ in slices: slice_.finalize()
    meep_utils.notify(model.simulation_name)

## Get the resonant field decay waveform inside the cavity
if meep.my_rank() == 0:
    ## Convert to polar notation and save the time-domain record
    x, y = np.array(x), np.array(y)
    meep_utils.savetxt(fname=model.simulation_name+"_timedomain.dat", X=zip(x, np.abs(y), meep_utils.get_phase(y)), fmt="%.6e",
            header=model.parameterstring + "#x-column time [s]\n#column ampli\n#column phase\n")

    with open("./last_simulation_name.dat", "w") as outfile: outfile.write(model.simulation_name) 

meep.all_wait()         # Wait until all file operations are finished
