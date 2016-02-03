#!/usr/bin/env python
#-*- coding: utf-8 -*-
""" An example python-meep simulation of a dielectric sphere scattering a broadband impulse, 
illustrating the use of the convenient functions provided by meep_utils.py 
(c) 2014 Filip Dominec, see http://f.dominec.eu/meep for more information 

The simulated models of structures are stored in the `metamaterial_models.py' module; see its code for the
description of each structure and command-line parameters accepted. Several of these parameters are not 
passed to the model; see the definition of `process_param()' in `meep_utils.py'.
"""

import time, sys, os
import numpy as np
from scipy.constants import c, epsilon_0, mu_0

import meep_utils, meep_materials, metamaterial_models
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep

# Model selection
model_param = meep_utils.process_param(sys.argv[1:])
model = metamaterial_models.models[model_param.get('model', 'default')](**model_param)

## Initialize volume, structure and the fields according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, pml_axes=meep.Z)
f = meep.fields(s)
# Define the Bloch-periodic boundaries (any transversal component of k-vector is allowed)
f.use_bloch(meep.X, getattr(model, 'Kx', 0) / (-2*np.pi)) 
f.use_bloch(meep.Y, getattr(model, 'Ky', 0) / (-2*np.pi))

# Add the field source (see meep_utils for an example of how an arbitrary source waveform is defined)
if not getattr(model, 'frequency', None):       ## (temporal source shape)
    #src_time_type = meep.band_src_time(model.src_freq/c, model.src_width/c, model.simtime*c/1.1)
    src_time_type = meep.gaussian_src_time(model.src_freq/c, model.src_width/c)
else:
    src_time_type = meep.continuous_src_time(getattr(model, 'frequency', None)/c)
srcvolume = meep.volume(                    ## (spatial source shape)
        meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness),
        meep.vec( model.size_x/2,  model.size_y/2, -model.size_z/2+model.pml_thickness))
f.add_volume_source(meep.Ex, src_time_type, srcvolume)

## Define monitor planes, and the field output for visualisation (controlled by keywords in the 'comment' parameter)
monitor_options = {'size_x':model.size_x, 'size_y':model.size_y, 'resolution':model.resolution, 'Kx':getattr(model, 'Kx', 0), 'Ky':getattr(model, 'Ky', 0)}
monitor1_Ex = meep_utils.AmplitudeMonitorPlane(f, comp=meep.Ex, z_position=model.monitor_z1, **monitor_options)
monitor1_Hy = meep_utils.AmplitudeMonitorPlane(f, comp=meep.Hy, z_position=model.monitor_z1, **monitor_options)
monitor2_Ex = meep_utils.AmplitudeMonitorPlane(f, comp=meep.Ex, z_position=model.monitor_z2, **monitor_options)
monitor2_Hy = meep_utils.AmplitudeMonitorPlane(f, comp=meep.Hy, z_position=model.monitor_z2, **monitor_options)

slices = []
if not "noepssnapshot" in model.comment:
    slices += [meep_utils.Slice(model=model, field=f, components=(meep.Dielectric), at_t=0, name='EPS')]
if "narrowfreq-snapshots" in model.comment:
    slices += [meep_utils.Slice(model=model, field=f, components=meep.Ex, at_y=0, at_t=np.inf,
            name=('At%.3eHz'%getattr(model, 'frequency', None)) if getattr(model, 'frequency', None) else '',
            outputpng=True, outputvtk=False)]
if "fieldevolution" in model.comment:
    slices += [meep_utils.Slice(model=model, field=f, components=(meep.Ex), at_x=0, name='FieldEvolution', 
        min_timestep=.1/model.src_freq, outputgif=True, outputhdf=True)]
if "snapshote" in model.comment:
    slices += [meep_utils.Slice(model=model, field=f, components=(meep.Ex, meep.Ey, meep.Ez), at_t=np.inf, name='SnapshotE')]

## Run the FDTD simulation or the frequency-domain solver
if not getattr(model, 'frequency', None):       ## time-domain computation
    f.step(); timer = meep_utils.Timer(simtime=model.simtime); meep.quiet(True) # use custom progress messages
    while (f.time()/c < model.simtime):     # timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
        for slice_ in slices: slice_.poll(f.time()/c)
    for slice_ in slices: slice_.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:                                       ## frequency-domain computation
    f.solve_cw(getattr(model, 'MaxTol',0.001), getattr(model, 'MaxIter', 5000), getattr(model, 'BiCGStab', 8)) 
    for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
    for slice_ in slices: slice_.finalize()
    meep_utils.notify(model.simulation_name)

## Get the reflection and transmission of the structure
if meep.my_rank() == 0:
    freq, s11, s12, columnheaderstring = meep_utils.get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, 
            frequency_domain=True if getattr(model, 'frequency', None) else False, 
            frequency=getattr(model, 'frequency', None),     ## procedure compatible with both FDTD and FDFD
            intf=getattr(model, 'interesting_frequencies', [0, model.src_freq+model.src_width]),  ## clip the frequency range for plotting
            pad_zeros=1.0,                                                                        ## speed-up FFT, and stabilize eff-param retrieval
            Kx=getattr(model, 'Kx', 0), Ky=getattr(model, 'Ky', 0),                                 ## enable oblique incidence (works only if monitors in vacuum)
            eps1=getattr(model, 'mon1eps', 1), eps2=getattr(model, 'mon2eps', 1))               ## enable monitors inside dielectrics

    meep_utils.savetxt(fname=model.simulation_name+".dat", fmt="%.6e",                            
            X=zip(freq, np.abs(s11), np.angle(s11), np.abs(s12), np.angle(s12)),                  ## Save 5 columns: freq, amplitude/phase for reflection/transmission
            header=model.parameterstring+columnheaderstring)     ## Export header

    with open("./last_simulation_name.dat", "w") as outfile: outfile.write(model.simulation_name) 

meep.all_wait()         # Wait until all file operations are finished
