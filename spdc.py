#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" run_sim.py - an example python-meep simulation of a dielectric sphere scattering a broadband impulse, 
illustrating the use of the convenient functions provided by meep_utils.py 
(c) 2014 Filip Dominec, see http://fzu.cz/~dominecf/meep/ for more information """


import numpy as np
import time, sys, os
import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep
c = 2.997e8

sim_param, model_param = meep_utils.process_param(sys.argv[1:])
class spdc_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=15e-12, resolution=3e-6, size_x=1350e-6, size_y=1350e-6, size_z=0,
            wgwidth=10e-6, wgheight=20e-6, monzd=180e-6):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "SPDC"    
        monzd=size_z

        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        substrate_z = size_x / 3
        self.pml_thickness = 10e-6
        self.monitor_z1, self.monitor_z2 = (-(monzd)/2, (monzd)/2)  
        self.simtime = simtime      # [s]
        self.srcFreq, self.srcWidth = 5000e9, 5000e9     # [Hz] (note: gaussian source ends at t=10/srcWidth)
        self.interesting_frequencies = (0e9, 2000e9)     # Which frequencies will be saved to disk
        self.Kx = 0; self.Ky = 0; self.padding=0
        self.size_x = size_x 
        self.size_y = size_y
        self.size_z = size_z

        ## Define materials
        self.materials   = [meep_materials.material_dielectric(eps=4., where = self.where_diel)]  
        #self.materials  += [meep_materials.material_dielectric(eps=4., where = self.where_substr)]  

        self.TestMaterials()
        f_c = c / np.pi/self.resolution/meep_utils.meep.use_Courant()
        meep_utils.plot_eps(self.materials, mark_freq=[f_c])

    # each material has one callback, used by all its polarizabilities (thus materials should never overlap)
    def where_diel(self, r):
        #curve1 =   self.size_x/np.pi * np.tanh(r.y()/self.size_y*3)
        #curve2 = - self.size_x/np.pi * np.tanh(r.y()/self.size_y*3)
        #if r.x() > curve1-self.wgwidth/2 and r.x() < curve1+self.wgwidth/2:
            #return self.return_value            # (do not change this line)
        #if r.x() > curve2-self.wgwidth/2 and r.x() < curve2+self.wgwidth/2 and r.y()>0:
            #return self.return_value            # (do not change this line)
        return 0
#}}}

# Model selection
model = spdc_model(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume and structure according to the model
#XXX vol = meep.vol2d(model.size_x, model.size_y, 1./model.resolution)
vol = meep.vol2d(model.size_x, model.size_y, 1./model.resolution)
vol.center_origin()
#s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes=meep.Z)
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes="All")

## Create fields with Bloch-periodic boundaries (any transversal component of k-vector is allowed, but may not radiate)
f = meep.fields(s)

## Add a source of the plane wave (see meep_utils for definition of arbitrary source shape)
if not sim_param['frequency_domain']:           ## Select the source dependence on time
    src_time_type = meep.band_src_time(model.srcFreq/c / 2 , model.srcWidth/c, model.simtime*c/1.1)
    #src_time_type = meep.gaussian_src_time(model.srcFreq/c, model.srcWidth/c)
else:
    src_time_type = meep.continuous_src_time(sim_param['frequency']/c)

# XXX srcvolume = meep.volume( 
        #meep.vec(-model.wgheight/2, -model.size_y/4-model.wgwidth/2, -model.size_z/2+model.pml_thickness),
        #meep.vec(+model.wgheight/2, -model.size_y/4+model.wgwidth/2, -model.size_z/2+model.pml_thickness))
srcvolume = meep.volume( 
        meep.vec(-model.size_x/2, -model.size_y/2+model.pml_thickness),
        meep.vec( model.size_x/2, -model.size_y/2+model.pml_thickness))
## Replace the f.add_volume_source(meep.Ex, srctype, srcvolume) line with following:
## Option for a custom source (e.g. exciting some waveguide mode)
class SrcAmplitudeFactor(meep.Callback): 
    ## The source amplitude is complex -> phase factor modifies its direction
    ## todo: implement in MEEP: we should define an AmplitudeVolume() object and reuse it for monitors later
    def __init__(self, Kx=0, Ky=0): meep.Callback.__init__(self)
    def complex_vec(self, vec):   ## Note: the 'vec' coordinates are _relative_ to the source center
        return (np.random.random()-.5) + 1j*(np.random.random()-.5)
af = SrcAmplitudeFactor(Kx=model.Kx, Ky=model.Ky)
meep.set_AMPL_Callback(af.__disown__())
f.add_volume_source(meep.Ez, src_time_type, srcvolume, meep.AMPL)

## Secondary (pump) source
src_time_type = meep.continuous_src_time(model.srcFreq/c)
f.add_volume_source(meep.Ez, src_time_type2, srcvolume)


## Define monitors planes and visualisation output
#monitor_options = {'size_x':model.size_x, 'size_y':model.size_y, 'Kx':model.Kx, 'Ky':model.Ky}
#monitor1_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z1, **monitor_options)
#monitor1_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z1, **monitor_options)
#monitor2_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z2, **monitor_options)
#monitor2_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z2, **monitor_options)

#XXX TODO 
slice_makers =  [meep_utils.Slice(model=model, field=f, components=(meep.Dielectric), at_t=0, name='EPS')]
slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ez, at_t=[0e-12, 100e-12], min_timestep=.025e-12, outputgif=True)]
slice_makers += [meep_utils.Slice(model=model, field=f, components=meep.Ez, at_t=2.5e-12)]

if not sim_param['frequency_domain']:       ## time-domain computation
    f.step()
    dt = (f.time()/c)
    meep_utils.lorentzian_unstable_check_new(model, dt)
    timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    while (f.time()/c < model.simtime):                               # timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        #for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
        for slice_maker in slice_makers: slice_maker.poll(f.time()/c)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    f.step()
    f.solve_cw(sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab']) 
    #for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name)

## Get the reflection and transmission of the structure
#if meep.my_rank() == 0:
    #freq, s11, s12 = meep_utils.get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, 
            #frequency_domain=sim_param['frequency_domain'], frequency=sim_param['frequency'], 
            #maxf=model.srcFreq+model.srcWidth, pad_zeros=1.0, Kx=model.Kx, Ky=model.Ky)
    #meep_utils.savetxt(freq=freq, s11=s11, s12=s12, model=model)
    #import effparam        # process effective parameters for metamaterials

with open("./last_simulation_name.txt", "w") as outfile: outfile.write(model.simulation_name) 
meep.all_wait()         # Wait until all file operations are finished
