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

sim_param, model_param = meep_utils.process_param(sys.argv[1:])
class HollowCyl_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=30e-9, resolution=5e-3, cellnumber=1, padding=9e-3, radius=33.774e-3, height=122.36e-3 ,  Kx=0, Ky=0): ## XXXheight=122.3642686e-3
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
    ## Convert to polar notation and save the time-domain record
    x, y = np.array(x), np.array(y)
    meep_utils.savetxt(fname=model.simulation_name+"_timedomain.dat", X=zip(x, np.abs(y), meep_utils.get_phase(y)), fmt="%.6e",
            header=model.parameterstring + meep_utils.sim_param_string(sim_param) + "#x-column time [s]\n#column ampli\n#column phase\n")

    with open("./last_simulation_name.dat", "w") as outfile: outfile.write(model.simulation_name) 

meep.all_wait()         # Wait until all file operations are finished
