#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" A model of a aperture-type near-field microscope for the terahertz range
(c) 2014 Filip Dominec, see http://f.dominec.eu/meep for more information """

import time, sys, os
import numpy as np
from scipy.constants import c, epsilon_0, mu_0

import meep_utils, meep_materials, metamaterial_models
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep

# Model selection
sim_param, model_param = meep_utils.process_param(sys.argv[1:])
class ApertureSphere_model(meep_utils.AbstractMeepModel): #{{{
    def __init__(self, comment="", simtime=100e-12, resolution=3e-6, Kx=0, Ky=0, 
            spacing=75e-6, monzd=50e-6,                              # lateral simulation size, and the z-length left for whole structure
            apertured=5e-6, apertureth=5e-6, gaasth=2e-6,           # metal aperture (square hole size, metal thickness, gallium arsenide layer thickness)
            radius=10e-6, epsloss=0.01, spherey=0e-6, spherez=-14e-6, # dielectric sphere (radius and dielectric losses)
            wireth=4e-6, wirey=14e-6, wirez=-14e-6                    # metallic wire along x  (diameter, lateral position)
            ):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "ApertureSphere"    
        self.register_locals(locals())          ## Remember the parameters

        ## Constants for the simulation
        self.pml_thickness = 20e-6
        self.monitor_z1, self.monitor_z2 = -monzd/2, self.apertureth+self.gaasth
        print  "self.monitor_z1, self.monitor_z2", self.monitor_z1, self.monitor_z2
        self.simtime = simtime      # [s]
        self.src_freq, self.src_width = 1000e9, 2000e9      # [Hz] (note: gaussian source ends at t=10/src_width)
        self.interesting_frequencies = (0e9, 5000e9)        # Which frequencies will be saved to disk

        self.size_x = spacing 
        self.size_y = spacing
        self.size_z = monzd*2 + 2*self.pml_thickness

        ## Define materials
        #self.materials = [meep_materials.material_dielectric(where = self.where_GaAs, eps=10)]  
        #self.materials += [meep_materials.material_TiO2_THz(where = self.where_TiO2)]  
        self.materials = []
        if self.radius > 0: self.materials += [meep_materials.material_dielectric(where = self.where_TiO2, eps=94., loss=epsloss)]  
        #if not 'NoGaAs' in comment:  self.materials += [meep_materials.material_Metal_THz(where = self.where_metal) ] %TODO  add GaAs
        self.materials += [meep_materials.material_Metal_THz(where = self.where_metal) ] ## TODO USE GOLD

        ## Test the validity of the model
        meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        self.test_materials()

    # each material has one callback, used by all its polarizabilities (thus materials should never overlap)
    def where_metal(self, r):
        if in_zslab(r, cz=self.apertureth/2, d=self.apertureth) and not (in_yslab(r, cy=0e-6, d=self.apertured) and in_xslab(r, cx=0e-6, d=self.apertured)): 
            return self.return_value        
        if in_zslab(r, cz=self.wirez, d=self.wireth) and in_yslab(r, cy=self.wirey, d=self.wireth): 
            return self.return_value     
        return 0

    def where_TiO2(self, r):
        if  in_sphere(r, cx=0, cy=self.spherey, cz=self.spherez, rad=self.radius):
            return self.return_value          
        return 0

    def where_GaAs(self, r):
        if in_zslab(r, cz=self.apertureth+(self.gaasth/2), d=self.gaasth):
            return self.return_value            
        return 0
#}}}

# Model selection
model = ApertureSphere_model(**model_param)
if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])

## Initialize volume, structure and the fields according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, sim_param=sim_param, pml_axes=meep.Z)
f = meep.fields(s)
f.use_bloch(meep.X, sim_param.get('Kx', 0) / (-2*np.pi)) # (any transversal component of k-vector is allowed)
f.use_bloch(meep.Y, sim_param.get('Ky',.0) / (-2*np.pi)) #  (TODO remove periodicity?)

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
monitor_options = {'size_x':model.size_x, 'size_y':model.size_y, 'Kx':0, 'Ky':0}
monitor1_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z1, **monitor_options)
monitor1_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z1, **monitor_options)
monitor_options = {'size_x':model.apertured, 'size_y':model.apertured, 'Kx':0, 'Ky':0} ## (specific for the apertured microscope)
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

    if not os.path.isfile('ref.dat'):   ## no reference yet, let us save one
        print "Saving the fields as a reference"
        meep_utils.savetxt(fname=model.simulation_name+".dat", fmt="%.6e",
                X=zip(freq, np.abs(s11), np.angle(s11), np.abs(s12), np.angle(s12)), 
                header=model.parameterstring + meep_utils.sim_param_string(sim_param) + headerstring)
    else:           ## save fields normalized to the reference
        print "Saving fields normalized to the reference (loaded from ref.dat)"
        (fref, s11refabs, s11refangle, s12refabs, s12refangle) = np.loadtxt('ref.dat', usecols=list(range(5)), unpack=True)
        meep_utils.savetxt(fname=model.simulation_name+"_NORMALIZED.dat", fmt="%.6e",
                X=zip(freq, np.abs(s11)/s11refabs, np.angle(s11)-s11refangle, np.abs(s12)/s12refabs, np.angle(s12)-s12refangle), 
                header=model.parameterstring + meep_utils.sim_param_string(sim_param) + headerstring)

    with open("./last_simulation_name.dat", "w") as outfile: outfile.write(model.simulation_name) 

meep.all_wait()         # Wait until all file operations are finished
