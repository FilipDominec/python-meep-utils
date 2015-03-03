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
from scipy.constants import pi, c
import meep_utils

import meep_mpi as meep
#import meep

meep.master_printf("=== Initialisation ===\n")
sim_param, model_param = meep_utils.process_param(sys.argv[1:])
print model_param

## Model selection
#from model_SphereWireNew import *;     model = SphereWire_model(**model_param) ## OK
#from model_SphereWireNew import *;     model = SphereFishnet_model(**model_param) ## OK
#from model_SphereWireNew import *;     model = SphereElliptic_model(**model_param)
#from model_SphereWireNew import *;     model = SimpleEllipsoid_model(**model_param)

#from model_CKEBars import *       
#model = CKEBars_model_test(**model_param)
#from model_PKCutSheet import *       
#model = PKCutSheet_model(**model_param)

from model_simple_structures import *       
#model = XCylWire_model(**model_param)
#model = YCylWire_model(**model_param)
#model = XRectWire_model(**model_param)
#model = XRectWireMet_model(**model_param)
#model = XCylWire_model_test(**model_param)
#model = dielbar_model(**model_param);      #  meep.use_averaging(True)
#model = PKCutSheet_model_test(**model_param)
#model = Fishnet_model(**model_param)
model = Wedge_model(**model_param)

#from model_SapphireBars import *       
#model = SapphireBars(**model_param)



if sim_param['frequency_domain']: model.simulation_name += ("_frequency=%.4e" % sim_param['frequency'])
meep.master_printf("Simulation name:\n\t%s\n" % model.simulation_name) ## TODO print parameters in a table

## Initialize volume
vol = meep.vol3d(model.size_x, model.size_y, model.size_z, 1./model.resolution)
volume_except_pml = meep.volume(
                meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness*0), 
                meep.vec(model.size_x/2,   model.size_y/2,  model.size_z/2-model.pml_thickness*0))
vol.center_origin()

## Define the Perfectly Matched Layers
perfectly_matched_layers = meep.pml(model.pml_thickness)          ## PML on both faces at Z axis

if not sim_param['frequency_domain']:
    meep.master_printf("== Time domain structure setup ==\n")
    ## Define each polarizability by redirecting the callback to the corresponding "where_material" function
    ## Define the frequency-independent epsilon for all materials (needed here, before defining s, or unstable)
    model.double_vec = model.eps; meep.set_EPS_Callback(model.__disown__())
    s = meep.structure(vol, meep.EPS, perfectly_matched_layers, meep.identity())

    ## Add all the materials
    model.build_polarizabilities(s)

    ## Add the source dependence
    #srctype = meep.band_src_time(model.srcFreq/c, model.srcWidth/c, model.simtime*c/1.1)
    srctype = meep.gaussian_src_time(model.srcFreq/c, model.srcWidth/c) ## , 0, 1000e-12    ?? 

else:
    meep.master_printf("== Frequency domain structure setup (for frequency of %g Hz) ==\n" % sim_param['frequency'])
    if (model.Kx!=0 or model.Ky!=0): print "Warning: frequency-domain solver may be broken for nonperpendicular incidence"
    ## Frequency-domain simulation does not support dispersive materials yet. We must define each material by 
    ## using the nondispersive permittivity and the nondispersive conductivity 
    ## (both calculated from polarizabilities at given frequency)

    ## Define the frequency-independent epsilon for all materials (has to be done _before_ defining s, or unstable)
    my_eps = meep_utils.MyHiFreqPermittivity(model, sim_param['frequency'])
    meep.set_EPS_Callback(my_eps.__disown__())
    s = meep.structure(vol, meep.EPS, perfectly_matched_layers, meep.identity()) 

    ## Create callback to set nondispersive conductivity (depends on material polarizabilities and frequency)
    mc = meep_utils.MyConductivity(model, sim_param['frequency'])
    meep.set_COND_Callback(mc.__disown__())
    s.set_conductivity(meep.E_stuff, meep.COND)  ## only "E_stuff" worked here for me
    srctype = meep.continuous_src_time(sim_param['frequency']/c)


## Create fields with Bloch-periodic boundaries (any nonzero transversal component of k-vector is possible)
f = meep.fields(s)
f.use_bloch(meep.X, -model.Kx/(2*np.pi))
f.use_bloch(meep.Y, -model.Ky/(2*np.pi))

## Add a source of a plane wave (with possibly oblique incidence)
## Todo implement in MEEP: we should define an AmplitudeVolume() object and reuse it for monitors later
srcvolume = meep.volume( 
        meep.vec(-model.size_x/2, -model.size_y/2, -model.size_z/2+model.pml_thickness),  ## TODO try from -inf to +inf
        meep.vec(model.size_x/2, model.size_y/2, -model.size_z/2+model.pml_thickness))
## TODO  move whole amplitude factor to meep_utils, exp(-1j*(a*x+b*y) - ((c*x)**2 + (d*y)**2))
class AmplitudeFactor(meep.Callback): 
    def __init__(self, Kx=0, Ky=0): 
        meep.Callback.__init__(self)
        (self.Kx, self.Ky) = Kx, Ky
    def complex_vec(self, vec):   ## Note: the 'vec' coordinates are _relative_ to the source center
        ## The source amplitude is complex and has the form of an oblique plane wave
        return np.exp(-1j*(self.Kx*vec.x() + self.Ky*vec.y()) - (vec.x()/.5e-3)**2 - (vec.y()/.5e-3)**2)
af = AmplitudeFactor(Kx=model.Kx, Ky=model.Ky)
meep.set_AMPL_Callback(af.__disown__())
f.add_volume_source(meep.Ex, srctype, srcvolume, meep.AMPL)
#f.add_volume_source(meep.Ex, srctype, srcvolume)

## Define monitors and visualisation output
monitor_options = {'size_x':model.size_x, 'size_y':model.size_y, 'Kx':model.Kx, 'Ky':model.Ky}
monitor1_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z1, **monitor_options)
monitor1_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z1, **monitor_options)
monitor2_Ex = meep_utils.AmplitudeMonitorPlane(comp=meep.Ex, z_position=model.monitor_z2, **monitor_options)
monitor2_Hy = meep_utils.AmplitudeMonitorPlane(comp=meep.Hy, z_position=model.monitor_z2, **monitor_options)
snapshot_maker = meep_utils.SnapshotMaker(snapshot_times=[model.simtime-float(X)/4/model.srcFreq for X in range(1)], 
        field=f, outputdir=model.simulation_name, volume=volume_except_pml)
#snapshot_maker = meep_utils.SnapshotMaker(snapshot_times=[], 
        #field=f, outputdir=model.simulation_name, volume=volume_except_pml)

#slice_makers = [meep_utils.SliceMaker(field=f, component=meep.Ex, timestep=.1e-12, normal="x", 
    #position=0., model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputHDF=True, outputVTK=True, outputGIF=False)]
#slice_makers = [meep_utils.SliceMaker(field=f, component=meep.Ex, timestep=.1e-12, normal="y", 
    #position=0., model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputHDF=True, outputVTK=True, outputGIF=False)]


pad = model.pml_thickness
slice_makers = [
        #meep_utils.SliceMaker(field=f, component=meep.Ex, timestep=.1e-12, normal="z", 
    #position=model.metalpos-2e-6, model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputHDF=True, outputVTK=True, outputGIF=False),

        #meep_utils.SliceMaker(field=f, component=meep.Ex, timestep=.1e-12, normal="z", 
    #position=15e-6, model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputHDF=True, outputVTK=True, outputGIF=False),
        #meep_utils.SliceMaker(field=f, component=meep.Ex, timestep=.1e-12, normal="z", 
    #position=20e-6, model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputHDF=True, outputVTK=True, outputGIF=False),
        #meep_utils.SliceMaker(field=f, component=meep.Ex, timestep=.1e-12, normal="z", 
    #position=25e-6, model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputHDF=True, outputVTK=True, outputGIF=False),

        #meep_utils.SliceMaker(field=f, component=meep.Ex, timestep=.1e-12, normal="x", 
    #position=0e-6, model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputHDF=True, outputVTK=True, outputGIF=False),

        ## 1D record - for the wedge numerical experiment
        meep_utils.SliceMaker(field=f, component=meep.Ex, timestep=.1e-12, 
                volume=meep.volume( 
                        meep.vec(0, -model.size_y/2+pad,  model.size_z/2-model.pml_thickness), 
                        meep.vec(0,  model.size_y/2-pad,  model.size_z/2-model.pml_thickness)),
                model=model, outputdir=model.simulation_name, pad=model.pml_thickness, outputHDF=True, outputVTK=True, outputGIF=True),

        ]
#slice_makers = []

meep.master_printf("=== Starting computation ===\n")

if not sim_param['frequency_domain']:
    tmptime = time.time()
    f.step()
    print 'setup took ', -tmptime + time.time() , 's'
    dt = (f.time()/c)
    meep_utils.lorentzian_unstable_check_new(model, dt)
    timer = meep_utils.Timer(simtime=model.simtime)
    #meep.quiet(True)
    count = 0
    while (f.time()/c < model.simtime):     ## timestepping cycle
        f.step()
        timer.print_progress(f.time()/c)
        for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
        for slice_maker in slice_makers: slice_maker.poll(f.time()/c)
        snapshot_maker.poll(f.time()/c)
        #print f.get_field(meep.Ex, meep.vec(0,0,0))
    meep.all_wait() ## FIXME needed?
    for slice_maker in slice_makers: slice_maker.finalize()
    meep_utils.notify(model.simulation_name, run_time=timer.get_time())
else:
    f.step()
    print sim_param['MaxIter']
    f.solve_cw(sim_param['MaxTol'], sim_param['MaxIter'], sim_param['BiCGStab']) 
    for monitor in (monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy): monitor.record(field=f)
    snapshot_maker.take_snapshot(0)
    meep_utils.notify(model.simulation_name)

with open("./last_simulation_name.txt", "w") as outfile: outfile.write(model.simulation_name) 

meep.master_printf("=== Processing recorded fields ===\n")
## Get the reflection and transmission of the structure
meep.master_printf("   getting s-params\n")
import time
if meep.my_rank() == 0:
    time1 = time.time()
    freq, s11, s12 = meep_utils.get_s_parameters(monitor1_Ex, monitor1_Hy, monitor2_Ex, monitor2_Hy, 
            frequency_domain=sim_param['frequency_domain'], 
            frequency=sim_param['frequency'], 
            maxf=model.srcFreq+model.srcWidth, 
            pad_zeros=1.0,
            Kx=model.Kx,
            Ky=model.Ky)
            #side_wavenumber=2*pi*modenumber*1/model.size_y)
    print "S-parameter retrieval (FFT etc.) took", time.time()-time1, "s" 
    #meep.master_printf("   saving\n")
    meep_utils.savetxt(freq=freq, s11=s11, s12=s12, model=model)
    import effparam
    #meep.master_printf("   done.\n")

print "All processes finishing", meep.my_rank()

meep.all_wait()         # Wait until all file operations are finished
