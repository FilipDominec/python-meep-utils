#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""
rtsim.py
This is a simulation of a broadband electromagnetic pulse propagating through a structure. 
The structure is loaded from a module. After the simulation, its s-parameters are calculated and saved.

About this script:
 * Written in 2012-2013 by Filip Dominec (dominecf at the server of fzu.cz)
 * Being distributed under the GPL license, this script is free as speech after five beers. 
 * You are encouraged to use and modify it as you need. Feel free to write me if needed.
 * Hereby I thank to the MEEP/python_meep authors and people of meep mailing list who helped me a lot.

Features and conventions: 
  * 3D simulation with Bloch-periodic walls, 
  * simulation outputs to r/t spectra, .gif animation and 3-D visualisation,
  * structure defined by a slow but versatile Python callback,
  * easy switching between time-domain and frequency-domain simulation, 
  * realistic dispersive material models,
  * calculation of complex amplitude parameters s11 (reflection), s21 (transmission)
  * all units in metric system (native time unit is then (1 m)/c = 3.33 ns),
  * meep remains in its own namespace ("import meep" instead of "from meep import *")
  
"""

import numpy as np
import time, sys, os
import numpy as np
from scipy.constants import c, epsilon_0, mu_0

import meep_utils, meep_materials
from meep_utils import in_sphere, in_xcyl, in_ycyl, in_zcyl, in_xslab, in_yslab, in_zslab, in_ellipsoid
import meep_mpi as meep
#import meep

class Wedge_model(meep_utils.AbstractMeepModel): #{{{
    """  Array of circular metallic wires along electric field
    FD 2013-07-13
    """
    def cell_centers(self):
        """ Helper function for stacked multilayered metamaterials """
        return np.arange(-self.monzd*(self.cells-1)/2, self.monzd*(self.cells-1)/2+1e-12, self.monzd)
        ## alternatively: add surrounding two cells to compare the propagation _inside_ a metamaterial
        #return np.arange(-self.monzd*(self.cells+1)/2, self.monzd*(self.cells+1)/2+1e-12, self.monzd)
    def __init__(self, comment="", simtime=100e-12, resolution=3e-6, cells=1, monzc=0e-6, Kx=0, Ky=0, padding=50e-6,
            radius=10e-6, yspacing=100e-6, zspacing=100e-6, monzd=200e-6, epsilon=100, **other_args):
        meep_utils.AbstractMeepModel.__init__(self)        ## Base class initialisation
        self.simulation_name = "Wedge"    ## 
        print(other_args)
        self.register_locals(locals(), other_args)          ## Remember the parameters

        ## Initialization of materials used
        if 'TiO2' in comment:
            self.materials = [meep_materials.material_TiO2_THz(where = self.where_wire)]
        elif 'STO' in comment:
            self.materials = [meep_materials.material_STO_THz(where = self.where_wire)]
        elif 'DielLossless' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.0)]
        elif 'DielLoLoss' in comment:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0.005)]
        else:
            self.materials = [meep_materials.material_dielectric(where = self.where_wire, eps=epsilon, loss=0)]


        ## Dimension constants for the simulation
        #self.size_x, self.size_y, self.size_z = xyspacing, xyspacing, 400e-6+cells*self.monzd
        self.pml_thickness = 50e-6
        self.size_x = resolution/1.8
        self.size_y = 2200e-6
        self.size_z = 2000e-6+cells*self.monzd + 2*self.padding + 2*self.pml_thickness

        ## constants for the simulation
        (self.monitor_z1, self.monitor_z2) = (-(monzd*cells)/2+monzc - padding, (monzd*cells)/2+monzc + padding)  
        self.simtime = simtime      # [s]
        self.src_width = 3000e9     
        self.src_freq = 4e9 + self.src_width/2 + (Ky**2+Kx**2)**.5 / (2*np.pi/3e8)  ## cutoff for oblique incidence
        self.interesting_frequencies = (0., 2000e9) 

        ## Test the validity of the model
        #meep_utils.plot_eps(self.materials, plot_conductivity=True, 
                #draw_instability_area=(self.f_c(), 3*meep.use_Courant()**2), mark_freq={self.f_c():'$f_c$'})
        try:
            self.test_materials() ## catches (most) errors in structure syntax, before they crash the callback
        except:
            print(sys.exc())

    def where_wire(self, r):
        y,z = r.y(), r.z()
        if z < self.size_z*(.3) and z-y/2>self.size_z*(-.2): 
            #return self.return_value
            yy,zz = np.arccos(np.cos(y/self.yspacing*np.pi*2))*self.yspacing/(np.pi*2), np.arccos(np.cos(z/self.zspacing*np.pi*2))*self.zspacing/(np.pi*2)
            if (yy**2+zz**2)**.5 < self.radius:
                return self.return_value
        return 0
#}}}


# Model selection
model_param = meep_utils.process_param(sys.argv[1:])
model = Wedge_model(**model_param)

## Initialize volume, structure and the fields according to the model
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
vol.center_origin()
s = meep_utils.init_structure(model=model, volume=vol, pml_axes=meep.Z)


## Create fields without Bloch-periodic boundaries
f = meep.fields(s)




# Add the field source (see meep_utils for an example of how an arbitrary source waveform is defined)
if not getattr(model, 'frequency', None):       ## (temporal source shape)
    #src_time_type = meep.band_src_time(model.src_freq/c, model.src_width/c, model.simtime*c/1.1)
    src_time_type = meep.gaussian_src_time(model.src_freq/c, model.src_width/c)
else:
    src_time_type = meep.continuous_src_time(getattr(model, 'frequency', None)/c)
srcvolume = meep.volume(                    ## (spatial source shape)
        meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness),
        meep.vec( model.size_x/2,  model.size_y/2, -model.size_z/2+model.pml_thickness))

#f.add_volume_source(meep.Ex, src_time_type, srcvolume)
## Replace the f.add_volume_source(meep.Ex, srctype, srcvolume) line with following:
## Option for a custom source (e.g. exciting some waveguide mode)
class SrcAmplitudeFactor(meep.Callback): 
    ## The source amplitude is complex -> phase factor modifies its direction
    ## todo: implement in MEEP: we should define an AmplitudeVolume() object and reuse it for monitors later
    def __init__(self, Kx=0, Ky=0): 
        meep.Callback.__init__(self)
        (self.Kx, self.Ky) = Kx, Ky
    def complex_vec(self, vec):   ## Note: the 'vec' coordinates are _relative_ to the source center
        # (oblique) Gaussian beam
        return np.exp(-1j*(self.Kx*vec.x() + self.Ky*vec.y()) - (vec.x()/.5e-3)**2 - (vec.y()/.5e-3)**2)
af = SrcAmplitudeFactor(Kx=getattr(model, 'Kx', 0), Ky=getattr(model, 'Ky', 0))
meep.set_AMPL_Callback(af.__disown__())
f.add_volume_source(meep.Ex, src_time_type, srcvolume, meep.AMPL)




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
        min_timestep=.1/model.src_freq, outputgif=True, outputvtk=True)]
if "snapshote" in model.comment:
    slices += [meep_utils.Slice(model=model, field=f, components=(meep.Ex, meep.Ey, meep.Ez), at_t=np.inf, name='SnapshotE')]

slices += [meep_utils.Slice(field=f, components=(meep.Ex), min_timestep=.1e-12, 
                volume=meep.volume( 
                        meep.vec(0, -model.size_y/2+model.pml_thickness, 0* model.size_z/2-model.pml_thickness), 
                        meep.vec(0,  model.size_y/2-model.pml_thickness, 0* model.size_z/2-model.pml_thickness)),
                model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputhdf=True, outputvtk=True, outputgif=True)]

        #pad = model.pml_thickness
        # 1D record - for the wedge numerical experiment
        #meep_utils.SliceMaker(field=f, component=meep.Ex, timestep=.1e-12, 
                #volume=meep.volume( 
                        #meep.vec(0, -model.size_y/2+pad,  model.size_z/2-model.pml_thickness), 
                        #meep.vec(0,  model.size_y/2-pad,  model.size_z/2-model.pml_thickness)),
                #model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputHDF=True, outputVTK=True, outputGIF=True),

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
    #t = monitor1_Ex.get_time()
    #Ex1, Hy1, Ex2, Hy2 = [mon.get_field_waveform() for mon in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy)]

    freq, s11, s12, columnheaderstring = meep_utils.get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, 
            frequency_domain=True if getattr(model, 'frequency', None) else False, 
            frequency=getattr(model, 'frequency', None),     ## procedure compatible with both FDTD and FDFD
            intf=getattr(model, 'interesting_frequencies', [0, model.src_freq+model.src_width]),  ## clip the frequency range for plotting
            pad_zeros=1.0,                                                                        ## speed-up FFT, and stabilize eff-param retrieval
            Kx=getattr(model, 'Kx', 0), Ky=getattr(model, 'Ky', 0),                                 ## enable oblique incidence (works only if monitors in vacuum)
            eps1=getattr(model, 'mon1eps', 1), eps2=getattr(model, 'mon2eps', 1))               ## enable monitors inside dielectrics

    print "EVERYTHING OK"
    meep_utils.savetxt(fname=model.simulation_name+".dat", fmt="%.6e",                            
            X=zip(freq, np.abs(s11), np.angle(s11), np.abs(s12), np.angle(s12)),                  ## Save 5 columns: freq, amplitude/phase for reflection/transmission
            header=model.parameterstring+columnheaderstring)     ## Export header

    with open("./last_simulation_name.dat", "w") as outfile: outfile.write(model.simulation_name) 

meep.all_wait()         # Wait until all file operations are finished
